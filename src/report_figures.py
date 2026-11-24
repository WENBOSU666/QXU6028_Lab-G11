# -*- coding: utf-8 -*-
"""
report_figures.py — regenerate the report's Figures 1–9 in the *same panel layout*
from the same data and model as src/l1_forward.py, src/l2_inverse.py and src/l3_map.py.

    python src/report_figures.py        (or through reproduce.py)

Writes output/report/fig01_report.png ... fig09_report.png
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "output" / "report"
OUT.mkdir(parents=True, exist_ok=True)

# palette (matches the report)
NAVY = "#1b2a4a"
TEAL = "#167b8b"
ORANGE = "#c0602a"
GREY = "#8c8c8c"
BAND = "#e8c9a0"
TEAL_FILL = "#8fc6bd"

kB = 8.617333e-5
NA = 3.816e15
E_CHARGE = 1.602176634e-19

# L3 settings
SIDE, REGION, SPEC, LAMBDA_REF = 10, 6.0, 3.0, 1.33e-3
PITCH = REGION / SIDE
FACTORS = (0.01, 0.1, 1.0, 10.0, 100.0)
N_NOISE, SEED = 500, 0


def style(ax, title=None, xlabel=None, ylabel=None):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=10)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(labelsize=9)


def save(fig, name):
    fig.savefig(OUT / name, dpi=230, bbox_inches="tight")
    plt.close(fig)
    print(f"[report figs] {name}")


# =========================================================== L1
def l1_data():
    import pandas as pd

    params = json.loads((DATA / "L1_params.json").read_text(encoding="utf-8-sig"))
    df = pd.read_csv(DATA / "L1_data.csv", encoding="utf-8-sig")
    EF = df["EF_eV"].to_numpy(float)
    n_band = df["n_band_cm2"].to_numpy(float)
    rho = df["rho_ohm_sq"].to_numpy(float)
    T, c = params["T"], params["c"]
    EF_op, rho_op = params["EF_op_eV"], params["rho_ohm_sq"]

    def count(ef, lower, upper=2.5):
        energy = np.linspace(lower, upper, 20001)
        f = 1.0 / (np.exp(np.clip((energy - ef) / (kB * T), -60, 60)) + 1.0)
        y = c * energy * f
        integ = np.trapezoid(y, energy) if hasattr(np, "trapezoid") else np.trapz(y, energy)
        return float(NA * integ)

    n_fermi = np.array([count(ef, ef) for ef in EF])
    C = n_band / n_fermi
    mu = 1.0 / (E_CHARGE * n_fermi * rho)
    return dict(EF=EF, n_band=n_band, n_fermi=n_fermi, C=C, mu=mu, rho=rho,
                T=T, c=c, EF_op=EF_op, rho_op=rho_op, count=count)


def fig01(d):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.1), layout="constrained")
    k = int(np.argmin(abs(d["EF"] - d["EF_op"])))
    ax = axes[0]
    ax.plot(d["EF"], d["n_band"] / 1e12, "o-", color=NAVY, ms=4, lw=1.6, label="Band edge")
    ax.plot(d["EF"], d["n_fermi"] / 1e12, "s-", color=TEAL, ms=4, lw=1.6, label="Fermi level")
    ax.plot([d["EF_op"]], [d["n_band"][k] / 1e12], "D", color=ORANGE, ms=8,
            label="Operating point")
    ax.plot([d["EF_op"]], [d["n_fermi"][k] / 1e12], "D", color=ORANGE, ms=8)
    ax.axvline(d["EF_op"], color="gray", ls=":", lw=0.9)
    style(ax, "(a) Carrier counts", "Fermi level (eV)", r"Density ($10^{12}$ cm$^{-2}$)")
    ax.legend(fontsize=8.5, frameon=False, loc="upper left")

    ax = axes[1]
    ax.plot(d["EF"], d["C"], "o-", color=TEAL, ms=4, lw=1.6)
    ax.axvline(d["EF_op"], color="gray", ls=":", lw=0.9)
    for target in (2.0, 4.0):
        ax.axhline(target, color="gray", ls=":", lw=0.9)
        ax.plot(0.07276 if target == 2 else 0.15735, target, "o", mfc="white",
                mec=TEAL, ms=7, mew=1.4)
    ax.plot([d["EF_op"]], [np.interp(d["EF_op"], d["EF"], d["C"])], "D", color=ORANGE, ms=8)
    style(ax, "(b) Correction across the sweep", "Fermi level (eV)", r"$n_{band}/n_F$")
    save(fig, "fig01_report.png")


def fig02(d):
    EF_op, T, c = d["EF_op"], d["T"], d["c"]
    E = np.linspace(0.0, 0.38, 4000)
    f = 1.0 / (np.exp(np.clip((E - EF_op) / (kB * T), -60, 60)) + 1.0)
    dens = NA * c * E * f / 1e12       # 10^12 cm^-2 eV^-1

    trap = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    below = trap(dens[E <= EF_op], E[E <= EF_op])
    above = trap(dens[E > EF_op], E[E > EF_op])
    pct_below = 100 * below / (below + above)

    fig, ax = plt.subplots(figsize=(11.0, 3.4), layout="constrained")
    ax.plot(E, dens, color=NAVY, lw=1.6)
    ax.fill_between(E, dens, where=(E <= EF_op), color=BAND, alpha=0.85,
                    label=f"Below the Fermi level ({pct_below:.1f}% of band count)")
    ax.fill_between(E, dens, where=(E >= EF_op), color=TEAL_FILL, alpha=0.85,
                    label=f"Fermi-referenced count ({100 - pct_below:.1f}%)")
    ax.axvline(EF_op, color=ORANGE, ls="--", lw=1.4)
    ax.annotate(r"$E_{F,op}$", (EF_op, dens.max() * 0.92), xytext=(6, 0),
                textcoords="offset points", color=ORANGE, fontsize=12)
    style(ax, None, "Energy above the Dirac point (eV)",
          "Occupied states\n($10^{12}$ cm$^{-2}$ eV$^{-1}$)")
    ax.legend(fontsize=8.5, frameon=False, loc="upper right")
    save(fig, "fig02_report.png")


# =========================================================== L2
def l2_data():
    meta = json.loads((DATA / "L2_meta.json").read_text(encoding="utf-8-sig"))
    raw = np.genfromtxt(DATA / "L2_data.csv", delimiter=",", names=True)
    T = float(meta["T_K"])
    kT = kB * T
    idx = np.unique(raw["measurement_index"]).astype(int)
    reps = []
    for i in idx:
        rows = raw[raw["measurement_index"] == i]
        rows = rows[np.argsort(rows["repeat"])]
        reps.append(np.column_stack([rows["R_H_raw"], rows["R_xx_raw"]]))
    reps = np.array(reps)
    means = reps.mean(axis=1)

    def normalise(pairs):
        return pairs / np.array([np.max(np.abs(pairs[:, 0])), np.max(pairs[:, 1])])

    energy = np.linspace(-2.5, 2.5, 4000)
    ef_grid = np.linspace(-8 * kT, 9 * kT, 400)
    gaps = np.r_[0.0, np.arange(0.5, 8.001, 0.25)]
    f = 1.0 / (np.exp(np.clip((energy[None, :] - ef_grid[:, None]) / kT, -60, 60)) + 1.0)
    w = np.full(4000, energy[1] - energy[0])
    w[[0, -1]] *= 0.5
    n_k = np.where(energy[None, :] >= ef_grid[:, None], f, 0) * w
    p_k = np.where(energy[None, :] <= ef_grid[:, None], 1 - f, 0) * w
    lib = []
    for gap in gaps:
        g = np.maximum(np.abs(energy) - gap * kT / 2, 0)
        n, p = n_k @ g, p_k @ g
        rh, rs = (p - n) / (p + n) ** 2, 1 / (p + n)
        lib.append(normalise(np.column_stack([rh, rs])))
    lib = np.array(lib)

    def fit(pairs):
        xy = normalise(pairs)
        d2 = ((lib[:, None, :, :] - xy[None, :, None, :]) ** 2).sum(axis=-1)
        nearest = d2.argmin(axis=2)
        rms = np.sqrt(d2.min(axis=2).mean(axis=1))
        best = 1 + np.argmin(rms[1:])
        chosen = best if (rms[best] < 0.85 * rms[0] and gaps[best] >= 1) else 0
        return chosen, ef_grid[nearest[chosen]], rms, best, nearest[chosen]

    chosen, traj, rms, best, nearest = fit(means)

    rng = np.random.default_rng(20260923)
    b_gap, b_ef = [], []
    for _ in range(200):
        pairs = reps[np.arange(len(idx)), rng.integers(0, 3, len(idx))]
        b, ef, _, _, _ = fit(pairs)
        b_gap.append(gaps[b])
        b_ef.append(ef)
    lo, hi = np.percentile(b_gap, [16, 84])
    ef_lo, ef_hi = np.percentile(b_ef, [16, 84], axis=0)

    return dict(meta=meta, kT=kT, T=T, idx=idx, reps=reps, means=means,
                normalise=normalise, energy=energy, ef_grid=ef_grid, gaps=gaps,
                lib=lib, rms=rms, chosen=chosen, best=best,
                nearest=nearest, traj=traj, lo=lo, hi=hi, ef_lo=ef_lo, ef_hi=ef_hi,
                model="hard_gap" if chosen else "cone")


def fig03(d):
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), layout="constrained")
    for ax, col, title, ylabel in [
            (axes[0], 0, "(a) Signed Hall signal", r"$R_H$ ($10^{-12}$ raw units)"),
            (axes[1], 1, "(b) Longitudinal signal", r"$R_{xx}$ ($10^{-12}$ raw units)")]:
        for rep in range(3):
            ax.plot(d["idx"], d["reps"][:, rep, col] * 1e12, "-", color=GREY, lw=0.8,
                    alpha=0.55, label="Raw repeat" if rep == 0 else None)
        ax.plot(d["idx"], d["means"][:, col] * 1e12, color=TEAL, lw=2.0, label="Repeat mean")
        ax.axvspan(24, 25, color=BAND, alpha=0.8)
        style(ax, title, "Measurement index", ylabel)
        if col == 0:
            ax.axhline(0, color="gray", ls="--", lw=0.9)
            ax.annotate("24–25", (25, 0), xytext=(34, 34), textcoords="offset points",
                        fontsize=9, color=NAVY,
                        arrowprops=dict(arrowstyle="-", color=GREY, lw=0.9))
        else:
            ax.legend(fontsize=8, frameon=False, loc="upper right")
    save(fig, "fig03_report.png")


def fig04(d):
    fig = plt.figure(figsize=(9.6, 5.6), layout="constrained")
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.92])

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(d["lib"][0, :, 0], d["lib"][0, :, 1], "--", color=GREY, lw=1.3, label="Cone")
    ax.plot(d["lib"][d["chosen"], :, 0], d["lib"][d["chosen"], :, 1], "-", color=TEAL,
            lw=1.7, label="Hard gap")
    xy = d["normalise"](d["means"])
    ax.plot(xy[:, 0], xy[:, 1], "o", mfc="none", mec=ORANGE, ms=5, mew=1.0,
            label="Paired means")
    style(ax, "(a) Hall-resistance fit", r"Normalised $R_H$", r"Normalised $R_{xx}$")
    ax.legend(fontsize=8, frameon=False, loc="center", bbox_to_anchor=(0.5, 0.45))

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(d["gaps"][1:], d["rms"][1:], "o-", color=TEAL, ms=3.5, lw=1.4)
    ax.plot(0, d["rms"][0], "s", color="#3b4a6b", ms=7, label="Gapless cone")
    ax.plot(d["gaps"][d["best"]], d["rms"][d["best"]], "D", color=ORANGE, ms=8)
    ax.annotate(f"{d['gaps'][d['best']]:g}", (d["gaps"][d["best"]], d["rms"][d["best"]]),
                xytext=(-26, 10), textcoords="offset points", color=ORANGE, fontsize=10)
    style(ax, "(b) Gap search", r"Trial gap ($k_B T$)", "RMS distance")
    ax.legend(fontsize=8, frameon=False, loc="upper right")

    ax = fig.add_subplot(gs[1, :])
    ax.fill_between(d["idx"], d["ef_lo"] * 1e3, d["ef_hi"] * 1e3, color=TEAL_FILL,
                    alpha=0.55, label="16th–84th percentile")
    ax.plot(d["idx"], d["traj"] * 1e3, "o-", color=NAVY, ms=3.5, lw=1.3, label="Fit to means")
    ax.axhline(0, color="gray", ls=":", lw=0.9)
    style(ax, "(c) Fermi level in measurement order", "Measurement index", r"$E_F$ (meV)")
    ax.legend(fontsize=8, frameon=False, loc="upper left", ncol=2)
    save(fig, "fig04_report.png")


def fig05_06(d):
    example = int(np.argmin(abs(d["gaps"] - 2)))
    xy = d["normalise"](d["means"])
    ef_meas = d["traj"] / d["kT"]
    best = d["gaps"][d["chosen"]]
    ylabels = [r"Normalised $R_H$", r"Normalised $R_s$"]
    for col, name, title in [
            (0, "fig05_report.png", "Hall coefficient versus Fermi level"),
            (1, "fig06_report.png", "Sheet resistance versus Fermi level")]:
        fig, ax = plt.subplots(figsize=(8.6, 4.6), layout="constrained")
        ax.plot(d["ef_grid"] / d["kT"], d["lib"][0, :, col], "--", color="#5b6b8c",
                lw=1.6, label=r"Cone, $E_g=0$")
        ax.plot(d["ef_grid"] / d["kT"], d["lib"][example, :, col], "-.", color=ORANGE,
                lw=1.6, label=r"Trial gap, $2.00\,k_BT$")
        ax.plot(d["ef_grid"] / d["kT"], d["lib"][d["chosen"], :, col], "-", color=TEAL,
                lw=1.9, label=rf"Best gap, ${best:.2f}\,k_BT$")
        ax.plot(ef_meas, xy[:, col], "o", mfc="none", mec="#2f4a5a", ms=4.2, mew=0.9,
                label=r"G11 means at fitted $E_F$")
        ax.axvline(0, color="gray", ls=":", lw=0.9)
        if col == 0:
            ax.axhline(0, color="gray", lw=0.8)
        style(ax, title, r"Fermi level, $E_F/k_BT$", ylabels[col])
        ax.legend(fontsize=9, frameon=False, loc="upper right")
        save(fig, name)


def fig07(d):
    xy = d["normalise"](d["means"])
    model = d["lib"][d["chosen"], d["nearest"]]
    res = xy - model
    rmse = np.sqrt((res ** 2).mean(axis=0))
    ef = d["traj"] / d["kT"]

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 5.4), layout="constrained")
    for ax, col, title in [(axes[0], 0, "(a) Hall coefficient"),
                           (axes[1], 1, "(b) Sheet resistance")]:
        c = TEAL if col == 0 else ORANGE
        ax.vlines(ef, 0, res[:, col], color=c, lw=0.9)
        ax.plot(ef, res[:, col], "o", color=c, ms=4)
        ax.axhline(0, color=NAVY, lw=0.9)
        ax.axvline(0, color="gray", ls=":", lw=0.9)
        ax.set_ylim(-0.04, 0.04)          # same fixed range as the report, both panels
        ax.set_yticks([-0.04, -0.02, 0.00, 0.02, 0.04])
        ax.text(0.99, 0.94, f"RMSE = {rmse[col]:.6f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=9.5)
        style(ax, title, None if col == 0 else r"Fitted $E_F/k_B T$", "Normalised residual")
        if col == 0:
            ax.tick_params(labelbottom=False)
    save(fig, "fig07_report.png")


# =========================================================== L3
def l3_setup():
    A = np.load(DATA / "L3_A.npy")
    D = np.loadtxt(DATA / "L3_d.csv", delimiter=",")
    n = A.shape[1]
    rows = []
    for r in range(SIDE):
        for c in range(SIDE):
            k = r * SIDE + c
            if c + 1 < SIDE:
                row = np.zeros(n); row[k], row[k + 1] = 1, -1; rows.append(row)
            if r + 1 < SIDE:
                row = np.zeros(n); row[k], row[k + SIDE] = 1, -1; rows.append(row)
    return A, D, np.asarray(rows)


def invert(A, d, lam, L):
    n = A.shape[1]
    m = np.linalg.solve(A.T @ A + lam * L.T @ L + 1e-12 * np.eye(n), A.T @ d)
    return m - m.mean()


def fig08(A, d, L):
    extent = [-3, 3, -3, 3]
    fig = plt.figure(figsize=(9.8, 3.5), layout="constrained")
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.15, 1.0, 1.15])

    grid = (invert(A, d, LAMBDA_REF, L) * 100).reshape(SIDE, SIDE)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(grid, extent=extent, origin="lower", cmap="RdBu_r", vmin=-3, vmax=3)
    for (r, c), mk in [(np.unravel_index(grid.argmax(), grid.shape), "o"),
                       (np.unravel_index(grid.argmin(), grid.shape), "s")]:
        ax.plot((c - 4.5) * PITCH, (r - 4.5) * PITCH, mk, ms=10, mfc="none", mec="black",
                mew=1.3)
    style(ax, "(a) Relative Hall conductivity", "x (probe pitches)", "y (probe pitches)")
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("Deviation from reference (%)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    sv = np.linalg.svd(A, compute_uv=False)
    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(np.arange(1, len(sv) + 1), sv, ".", ms=2.2, color="#1f3b73")
    style(ax, "(b) Forward matrix spectrum", "Singular value index", "Singular value")

    res, rough, labels = [], [], []
    for f in FACTORS:
        mf = invert(A, d, LAMBDA_REF * f, L)
        res.append(np.linalg.norm(A @ mf - d))
        rough.append(np.linalg.norm(L @ mf))
        labels.append(f"{f:g}×")
    ax = fig.add_subplot(gs[0, 2])
    ax.loglog(res, rough, "o-", ms=4, color=TEAL)
    for r, g, lab, off in zip(res, rough, labels,
                              [(4, 6), (4, -10), (6, 6), (6, 4), (-18, 6)]):
        ax.annotate(lab, (r, g), xytext=off, textcoords="offset points", fontsize=7)
    style(ax, "(c) Changing the penalty", r"Residual norm $\|Am-d\|_2$",
          r"Roughness $\|Lm\|_2$")
    save(fig, "fig08_report.png")


def fig09(A, D, L):
    extent = [-3, 3, -3, 3]
    fig = plt.figure(figsize=(9.8, 3.2), layout="constrained")
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 1.0, 1.35])
    for col, (r0, c0), tag in [(0, (5, 5), "(a) Central impulse"),
                               (1, (2, 7), "(b) Off-centre impulse")]:
        e = np.zeros(A.shape[1]); e[r0 * SIDE + c0] = 1.0
        psf = invert(A, A @ e, LAMBDA_REF, L).reshape(SIDE, SIDE)
        ax = fig.add_subplot(gs[0, col])
        ax.imshow(psf, extent=extent, origin="lower", cmap="viridis")
        ax.plot((c0 - 4.5) * PITCH, (r0 - 4.5) * PITCH, "x", color="white", ms=8, mew=1.4)
        style(ax, tag, "x (pitches)", "y (pitches)" if col == 0 else None)

    sigma = D.std(axis=0, ddof=1)
    rng = np.random.default_rng(SEED)
    peaks = [np.abs(invert(A, rng.normal(0.0, sigma), LAMBDA_REF, L)).max() * 1e6
             for _ in range(N_NOISE)]
    det = float(np.percentile(peaks, 95))
    ax = fig.add_subplot(gs[0, 2])
    ax.hist(peaks, bins=26, color="#7fb3a9", edgecolor="white")
    ax.axvline(det, color="#c0392b", ls="--", lw=1.2)
    ax.annotate(f"95th percentile\n{det:.2f}", (det, ax.get_ylim()[1] * 0.72), xytext=(6, 0),
                textcoords="offset points", fontsize=7.5)
    style(ax, "(c) Noise-only peak", r"Peak $|m|$ ($10^{-6}$)", "Count")
    save(fig, "fig09_report.png")
    print(f"[report figs] detection limit {det:.2f}e-6 (report 2.74)")


def main():
    d1 = l1_data()
    fig01(d1)
    fig02(d1)

    d2 = l2_data()
    print(f"[report figs] L2: model={d2['model']} Eg={d2['gaps'][d2['chosen']]:g} kT "
          f"rms {d2['rms'][0]:.6f} -> {d2['rms'][d2['best']]:.6f}")
    fig03(d2)
    fig04(d2)
    fig05_06(d2)
    fig07(d2)

    A, D, L = l3_setup()
    fig08(A, D.mean(axis=0), L)
    fig09(A, D, L)
    print(f"[report figs] all written to {OUT}")


if __name__ == "__main__":
    main()
