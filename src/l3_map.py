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
LAMBDA_REF = 1.33e-3       # reference smoothness strength
SWEEP_FACTORS = (0.01, 0.1, 1.0, 10.0, 100.0)
N_NOISE = 500
SEED = 0


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


class Inverter:
    def __init__(self, A: np.ndarray, lam: float):
        self.A = A
        self.n = A.shape[1]
        self.side = int(round(np.sqrt(self.n)))
        self.L = difference_operator(self.n, self.side)
        self.lhs = A.T @ A + lam * (self.L.T @ self.L) + 1e-12 * np.eye(self.n)

    def __call__(self, rhs: np.ndarray) -> np.ndarray:
        """Solve the regularised normal equations and project onto zero mean."""
        m = np.linalg.solve(self.lhs, rhs)
        return m - m.mean()


def build(lam: float):
    A = np.load(DATA / "L3_A.npy")
    D = np.loadtxt(DATA / "L3_d.csv", delimiter=",")      # (3, 128)
    meta = json.loads((DATA / "L3_meta.json").read_text(encoding="utf-8-sig"))
    side = int(meta.get("npix", N_PIX))
    pitch = REGION_IN_PITCH / side
    inv = Inverter(A, lam)
    d = D.mean(axis=0)
    sigma = float(D.std(axis=0, ddof=1).mean())
    return A, D, meta, side, pitch, inv, d, sigma


def psf_width(inv: Inverter, side: int, row: int, col: int) -> float:
    """Width of the reconstructed point-spread function, in pixels (2 sigma)."""
    e = np.zeros(inv.n)
    e[row * side + col] = 1.0
    m = inv(inv.A.T @ (inv.A @ e)).reshape(side, side)
    lobe = np.clip(m, 0, None)
    total = lobe.sum()
    if total <= 0:
        return float("nan")
    yy, xx = np.mgrid[0:side, 0:side]
    cy = (lobe * yy).sum() / total
    cx = (lobe * xx).sum() / total
    var = (lobe * ((yy - cy) ** 2 + (xx - cx) ** 2)).sum() / total
    return float(2.0 * np.sqrt(var))


def main():
    A, D, meta, side, pitch, inv, d, sigma = build(LAMBDA_REF)

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
        inv_f = Inverter(A, LAMBDA_REF * f)
        g = (inv_f(A.T @ d) * 100.0)
        sweep.append({"factor": f, "max_pct": float(g.max()),
                      "max_row": int(np.unravel_index(np.argmax(g), grid.shape)[0]),
                      "max_col": int(np.unravel_index(np.argmax(g), grid.shape)[1])})

    # 3. spatial resolution from the point-spread function
    w_centre = psf_width(inv, side, side // 2, side // 2)
    w_off = psf_width(inv, side, 2, side - 3)
    resolution_pitch = max(w_centre, w_off) * pitch

    # 4. detection limit: 500 noise-only reconstructions, 95th percentile
    rng = np.random.default_rng(SEED)
    peaks = []
    for _ in range(N_NOISE):
        nv = rng.normal(0.0, sigma, size=A.shape[0])
        peaks.append(np.abs(inv(A.T @ nv) * 100.0).max())
    detection_pct = float(np.percentile(peaks, 95))

    # 5. decision
    margin = float(np.abs(grid).max() + detection_pct)
    passed = bool(margin < SPEC_PERCENT)

    result = {
        "map": [[round(float(v), 6) for v in row] for row in grid],
        "features": [
            {"kind": "maximum", "row": int(r_max), "col": int(c_max),
             "x_pitch": round(x_max, 4), "y_pitch": round(y_max, 4),
             "deviation_pct": round(float(grid.max()), 4)},
            {"kind": "minimum", "row": int(r_min), "col": int(c_min),
             "x_pitch": round(x_min, 4), "y_pitch": round(y_min, 4),
             "deviation_pct": round(float(grid.min()), 4)},
        ],
        "resolution_probe_pitch": round(float(resolution_pitch), 3),
        "detection_limit_fraction": round(detection_pct / 100.0, 9),
        "pass": passed,
        "confidence": (f"95%: detection limit {detection_pct:.5f}% from {N_NOISE} noise-only "
                       f"reconstructions; margin {margin:.4f}% against {SPEC_PERCENT:.1f}%"),
        "diagnostics": {
            "lambda_ref": LAMBDA_REF,
            "pitch_probe_pitch": pitch,
            "sigma_per_reading": sigma,
            "margin_pct": margin,
            "resolution_centre_pitch": round(w_centre * pitch, 3),
            "resolution_offcentre_pitch": round(w_off * pitch, 3),
            "sweep": sweep,
        },
    }
    (OUT / "L3_results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("features", "resolution_probe_pitch",
                                             "detection_limit_fraction", "pass", "confidence")},
                     indent=2, ensure_ascii=False))
    return result


if __name__ == "__main__":
    main()
