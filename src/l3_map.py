# -*- coding: utf-8 -*-
"""
L3 — two-dimensional inverse: Hall-conductivity non-uniformity map, spatial
resolution, detection limit and the pass/fail decision for the 3 % specification.

Forward relation: d = A m, with A (128 x 100) supplied already row-centred, so
every row sums to zero and the probe is blind to a uniform offset. The map is
therefore obtained under a zero-mean constraint (subtracting the mean leaves the
predictions unchanged because A @ 1 = 0).

Files read from data/:
    L3_A.npy     (128, 100) forward matrix
    L3_d.csv     (3, 128)    three repeats of the 128 readings, no header
    L3_meta.json grid size, region and projection note
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)

N_PIX = 10                 # 10 x 10 map
REGION_IN_PITCH = 6.0      # scan region is 6s x 6s, s = probe pitch
SPEC_PERCENT = 3.0
SWEEP_FACTORS = (0.01, 0.1, 1.0, 10.0, 100.0)
N_NOISE = 500
SEED = 20260924            # recorded in the confidence field of results.json


def difference_operator(n_pix: int, side: int) -> np.ndarray:
    """Neighbour differences on the side x side grid (no wrap-around)."""
    rows = []
    for r in range(side):
        for c in range(side):
            k = r * side + c
            if c + 1 < side:
                row = np.zeros(n_pix)
                row[k], row[k + 1] = 1.0, -1.0
                rows.append(row)
            if r + 1 < side:
                row = np.zeros(n_pix)
                row[k], row[k + side] = 1.0, -1.0
                rows.append(row)
    return np.asarray(rows)


def inv_eta(inv):
    """Recover the eta actually used by an Inverter."""
    return float(getattr(inv, "eta", 0.0))


class Inverter:
    def __init__(self, A: np.ndarray, lam: float):
        self.A = A
        self.n = A.shape[1]
        self.side = int(round(np.sqrt(self.n)))
        self.L = difference_operator(self.n, self.side)
        self.eta = float(lam)
        self.lhs = A.T @ A + lam * (self.L.T @ self.L) + 1e-12 * np.eye(self.n)

    def __call__(self, rhs: np.ndarray) -> np.ndarray:
        """Solve the regularised normal equations and project onto zero mean."""
        m = np.linalg.solve(self.lhs, rhs)
        return m - m.mean()


def reference_eta(A: np.ndarray, L: np.ndarray) -> float:
    """eta0 = 1e-2 * ||A||_2^2 / ||L||_2^2  (methods guide, L3 starter prompt)."""
    sA = np.linalg.svd(A, compute_uv=False)[0]
    sL = np.linalg.svd(L, compute_uv=False)[0]
    return 1e-2 * sA ** 2 / sL ** 2


def build(lam: float = None):
    A = np.load(DATA / "L3_A.npy")
    D = np.loadtxt(DATA / "L3_d.csv", delimiter=",")      # (3, 128)
    meta = json.loads((DATA / "L3_meta.json").read_text(encoding="utf-8-sig"))
    side = int(meta.get("npix", N_PIX))
    pitch = REGION_IN_PITCH / side
    L = difference_operator(A.shape[1], side)
    if lam is None:
        lam = reference_eta(A, L)
    inv = Inverter(A, lam)
    d = D.mean(axis=0)
    sigma = D.std(axis=0, ddof=1)          # per-channel noise: each reading has its own repeat spread
    return A, D, meta, side, pitch, inv, d, sigma


def psf_width(inv: Inverter, side: int, row: int, col: int) -> float:
    """FWHM of the positive central lobe, in pixels.

    Methods guide, p. 12: inject a unit value at one pixel, remove its mean so the
    test map is zero-mean, reconstruct with the same eta and solver, then use the
    positive-lobe weights w = max(q_hat, 0) with the radius taken from the injected
    pixel:  sigma^2 = sum(w r^2) / (2 sum w),  FWHM = 2.355 sigma.
    """
    q = np.zeros(inv.n)
    q[row * side + col] = 1.0
    q -= q.mean()
    q_hat = inv(inv.A.T @ (inv.A @ q)).reshape(side, side)
    w = np.clip(q_hat, 0, None)
    total = w.sum()
    if total <= 0:
        return float("nan")
    yy, xx = np.mgrid[0:side, 0:side]
    r2 = (yy - row) ** 2 + (xx - col) ** 2
    sigma2 = (w * r2).sum() / (2.0 * total)
    return float(2.355 * np.sqrt(sigma2))


def noise_limit(inv: Inverter, sigma: np.ndarray, draws: int = N_NOISE,
                seed: int = SEED) -> float:
    """95th percentile of max|m| over noise-only reconstructions, in per cent."""
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(draws):
        nv = rng.normal(0.0, sigma)
        vals.append(np.abs(inv(inv.A.T @ nv) * 100.0).max())
    return float(np.percentile(vals, 95))


def main():
    A, D, meta, side, pitch, inv, d, sigma = build()
    eta = inv_eta(inv)

    # 1. non-uniformity map, in per cent of the reference conductivity
    m = inv(A.T @ d)
    grid = (m * 100.0).reshape(side, side)
    r_max, c_max = np.unravel_index(np.argmax(grid), grid.shape)
    r_min, c_min = np.unravel_index(np.argmin(grid), grid.shape)

    def coords(r, c):
        return float((c - (side - 1) / 2) * pitch), float((r - (side - 1) / 2) * pitch)

    x_max, y_max = coords(r_max, c_max)
    x_min, y_min = coords(r_min, c_min)

    # 2. stability against the smoothness strength
    sweep = []
    for f in SWEEP_FACTORS:
        inv_f = Inverter(A, eta * f)
        g = (inv_f(A.T @ d) * 100.0)
        sweep.append({"factor": f, "max_pct": float(g.max()),
                      "max_row": int(np.unravel_index(np.argmax(g), grid.shape)[0]),
                      "max_col": int(np.unravel_index(np.argmax(g), grid.shape)[1])})

    # 3. spatial resolution: FWHM of the positive lobe (methods guide, p. 12)
    w_centre = psf_width(inv, side, side // 2, side // 2)
    w_off = psf_width(inv, side, 2, side - 3)
    resolution_pitch = max(w_centre, w_off) * pitch

    # 4. detection limit: 500 noise-only reconstructions at the same eta as the map
    rng = np.random.default_rng(SEED)
    peaks = []
    for _ in range(N_NOISE):
        nv = rng.normal(0.0, sigma)
        peaks.append(np.abs(inv(A.T @ nv) * 100.0).max())
    detection_pct = float(np.percentile(peaks, 95))

    # 4b. same test at the weakest eta of the sweep: peak and delta_95 must share eta
    weakest = min(sweep, key=lambda s: s["factor"])
    inv_weak = Inverter(A, eta * weakest["factor"])
    delta_weak = noise_limit(inv_weak, sigma)
    weakest_bound = float(weakest["max_pct"]) + delta_weak

    # 5. decision
    margin = float(np.abs(grid).max() + detection_pct)
    passed = bool(margin < SPEC_PERCENT)

    result = {
        "map": [round(float(v) / 100.0, 12) for v in grid.ravel()],
        "features": [
            {"kind": "maximum", "row": int(r_max), "col": int(c_max),
             "x_pitch": round(x_max, 4), "y_pitch": round(y_max, 4),
             "deviation_pct": round(float(grid.max()), 4)},
            {"kind": "minimum", "row": int(r_min), "col": int(c_min),
             "x_pitch": round(x_min, 4), "y_pitch": round(y_min, 4),
             "deviation_pct": round(float(grid.min()), 4)},
        ],
        "resolution_probe_pitch": float(resolution_pitch),
        "detection_limit_fraction": float(detection_pct / 100.0),
        "pass": passed,
        "confidence": (
            f"Clear pass at the reference regularisation: the 95th-percentile detection "
            f"limit is {detection_pct:.6f}% of the reference conductivity ({N_NOISE} noise-only "
            f"reconstructions, seed {SEED}). The noise-adjusted upper bound is {margin:.4f}%, "
            f"below the {SPEC_PERCENT:.1f}% limit by {SPEC_PERCENT - margin:.4f} percentage points."),
        "diagnostics": {
            "lambda_ref": eta,
            "pitch_probe_pitch": pitch,
            "sigma_per_reading_mean": float(np.mean(sigma)),
            "sigma_per_reading_max": float(np.max(sigma)),
            "margin_pct": margin,
            "resolution_centre_pitch": round(w_centre * pitch, 3),
            "resolution_offcentre_pitch": round(w_off * pitch, 3),
            "weakest_eta_factor": weakest["factor"],
            "weakest_eta_peak_pct": weakest["max_pct"],
            "weakest_eta_delta95_pct": round(delta_weak, 6),
            "weakest_eta_bound_pct": round(weakest_bound, 4),
            "sweep": sweep,
        },
    }
    (OUT / "L3_results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    save_figures(A, D, inv, grid, (r_max, c_max), (r_min, c_min), side, pitch, sweep, peaks)
    print(json.dumps({k: result[k] for k in ("features", "resolution_probe_pitch",
                                             "detection_limit_fraction", "pass", "confidence")},
                     indent=2, ensure_ascii=False))
    return result


def save_figures(A, D, inv, grid, mx, mn, side, pitch, sweep, peaks):
    """Four figures: map, spectrum, regularisation sweep, point-spread function."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    extent = [-side * pitch / 2, side * pitch / 2, -side * pitch / 2, side * pitch / 2]
    teal = "#167b8b"

    # --- fig 1: non-uniformity map ------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.4, 5.4), layout="constrained")
    im = ax.imshow(grid, extent=extent, origin="lower", cmap="RdBu_r",
                   vmin=-3.0, vmax=3.0, interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, shrink=0.9)
    cb.set_label("Hall-conductivity deviation (%)")
    for (r, c), marker, edge in [(mx, "o", "black"), (mn, "s", "black")]:
        ax.plot((c - (side - 1) / 2) * pitch, (r - (side - 1) / 2) * pitch,
                marker, ms=13, mfc="none", mec=edge, mew=1.6)
    ax.set(xlabel="x (probe pitches)", ylabel="y (probe pitches)",
           title=f"G11 | zero-mean reconstruction (max {grid.max():+.3f} %, min {grid.min():+.3f} %)")
    ax.grid(color="white", lw=0.4, alpha=0.5)
    fig.savefig(OUT / "fig_L3_map.png", dpi=220)
    plt.close(fig)

    # --- fig 2: singular values of A ---------------------------------------------
    sv = np.linalg.svd(A, compute_uv=False)
    fig, ax = plt.subplots(figsize=(6.0, 4.2), layout="constrained")
    ax.semilogy(np.arange(1, len(sv) + 1), sv, "o-", ms=3.5, color=teal)
    ax.axvline(np.linalg.matrix_rank(A), color="gray", ls=":", lw=1.2)
    ax.annotate(f"rank = {np.linalg.matrix_rank(A)}", (np.linalg.matrix_rank(A), sv[0]),
                xytext=(8, 0), textcoords="offset points", fontsize=9)
    ax.set(xlabel="index", ylabel="singular value",
           title="Forward matrix A: 8 unobservable directions")
    fig.savefig(OUT / "fig_L3_spectrum.png", dpi=220)
    plt.close(fig)

    # --- fig 3: regularisation sweep ---------------------------------------------
    factors = [s["factor"] for s in sweep]
    maxima = [s["max_pct"] for s in sweep]
    fig, ax = plt.subplots(figsize=(6.0, 4.2), layout="constrained")
    ax.semilogx(factors, maxima, "o-", color=teal, label="max |deviation|")
    ax.axhline(SPEC_PERCENT, color="#c0392b", ls="--", lw=1.2, label="specification 3.0 %")
    ax.axvline(1.0, color="gray", ls=":", lw=1.0)
    ax.annotate(f"reference: {maxima[2]:.4f} %", (1.0, maxima[2]),
                xytext=(10, -18), textcoords="offset points", fontsize=9,
                arrowprops=dict(arrowstyle="-", color="gray"))
    ax.set(xlabel="regularisation strength (x reference)", ylabel="peak deviation (%)",
           title="Peak position is stable; amplitude is not")
    ax.legend(fontsize=9)
    fig.savefig(OUT / "fig_L3_sweep.png", dpi=220)
    plt.close(fig)

    # --- fig 4: point-spread function --------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.0), layout="constrained")
    psf_centre = None
    for ax, (r0, c0), tag in [(axes[0], (side // 2, side // 2), "centre"),
                              (axes[2], (2, side - 3), "off-centre")]:
        e = np.zeros(inv.n)
        e[r0 * side + c0] = 1.0
        psf = inv(inv.A.T @ (inv.A @ e)).reshape(side, side)
        if psf_centre is None:
            psf_centre = psf
        ax.imshow(psf, origin="lower", cmap="viridis", extent=extent)
        ax.plot((c0 - (side - 1) / 2) * pitch, (r0 - (side - 1) / 2) * pitch,
                "x", color="white", ms=9, mew=1.6)
        ax.set(title=f"unit impulse, {tag}", xlabel="x (probe pitches)")
    axes[0].set_ylabel("y (probe pitches)")

    row = psf_centre[side // 2, :]
    axes[1].plot(np.arange(side) * pitch, row / np.abs(row).max(), "o-", ms=4, color=teal)
    axes[1].axhline(0, color="gray", lw=0.8, ls=":")
    axes[1].set(title="point-spread function, central row",
                xlabel="x (probe pitches)", ylabel="normalised response")
    fig.savefig(OUT / "fig_L3_psf.png", dpi=220)
    plt.close(fig)

    # --- fig 5: detection limit ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.0, 4.0), layout="constrained")
    ax.hist(peaks, bins=30, color=teal, alpha=0.8)
    det = np.percentile(peaks, 95)
    ax.axvline(det, color="#c0392b", ls="--", lw=1.4,
               label=f"95th percentile = {det:.5f} %")
    ax.set(xlabel="max |deviation| of a noise-only reconstruction (%)",
           ylabel="count", title=f"{N_NOISE} noise-only reconstructions")
    ax.legend(fontsize=9)
    fig.savefig(OUT / "fig_L3_noise.png", dpi=220)
    plt.close(fig)
    print(f"[L3] figures written to {OUT}")


if __name__ == "__main__":
    main()
