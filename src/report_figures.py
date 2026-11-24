# -*- coding: utf-8 -*-
"""
report_figures.py — regenerate the figures *in the layout used in the report*,
from the same pipeline as src/l1_forward.py, src/l2_inverse.py and src/l3_map.py.

Writes to output/report/:
    fig08_report.png   (a) relative Hall conductivity  (b) forward matrix spectrum
                       (c) changing the penalty: roughness vs residual
    fig09_report.png   (a) central impulse  (b) off-centre impulse  (c) noise-only peak
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

N_PIX, SIDE = 100, 10
REGION = 6.0
PITCH = REGION / SIDE
SPEC = 3.0
LAMBDA_REF = 1.33e-3
FACTORS = (0.01, 0.1, 1.0, 10.0, 100.0)
N_NOISE = 500
SEED = 0


def load():
    A = np.load(DATA / "L3_A.npy")
    D = np.loadtxt(DATA / "L3_d.csv", delimiter=",")
    meta = json.loads((DATA / "L3_meta.json").read_text(encoding="utf-8-sig"))
    return A, D, meta


def diff_op(n, side):
    rows = []
    for r in range(side):
        for c in range(side):
            k = r * side + c
            if c + 1 < side:
                row = np.zeros(n); row[k], row[k + 1] = 1, -1; rows.append(row)
            if r + 1 < side:
                row = np.zeros(n); row[k], row[k + side] = 1, -1; rows.append(row)
    return np.asarray(rows)


def invert(A, d, lam, L):
    n = A.shape[1]
    m = np.linalg.solve(A.T @ A + lam * L.T @ L + 1e-12 * np.eye(n), A.T @ d)
    return m - m.mean()


def fig08(A, d, L):
    extent = [-3, 3, -3, 3]
    fig = plt.figure(figsize=(9.8, 3.5), layout="constrained")
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.15, 1.0, 1.15])

    m = invert(A, d, LAMBDA_REF, L)
    grid = (m * 100).reshape(SIDE, SIDE)
    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(grid, extent=extent, origin="lower", cmap="RdBu_r", vmin=-3, vmax=3)
    for (r, c), mk in [(np.unravel_index(grid.argmax(), grid.shape), "o"),
                       (np.unravel_index(grid.argmin(), grid.shape), "s")]:
        ax.plot((c - 4.5) * PITCH, (r - 4.5) * PITCH, mk, ms=10, mfc="none",
                mec="black", mew=1.3)
    ax.set(title="(a) Relative Hall conductivity", xlabel="x (probe pitches)",
           ylabel="y (probe pitches)")
    ax.title.set_fontsize(9)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("Deviation from reference (%)", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    sv = np.linalg.svd(A, compute_uv=False)
    ax = fig.add_subplot(gs[0, 1])
    ax.semilogy(np.arange(1, len(sv) + 1), sv, ".", ms=2.2, color="#1f3b73")
    ax.set(title="(b) Forward matrix spectrum", xlabel="Singular value index",
           ylabel="Singular value")

    res, rough, labels = [], [], []
    for f in FACTORS:
        mf = invert(A, d, LAMBDA_REF * f, L)
        res.append(np.linalg.norm(A @ mf - d))
        rough.append(np.linalg.norm(L @ mf))
        labels.append(f"{f:g}×" if f >= 1 else f"{f:g}×")
    ax = fig.add_subplot(gs[0, 2])
    ax.loglog(res, rough, "o-", ms=4, color="#167b8b")
    offsets = [(4, 6), (4, -10), (6, 6), (6, 4), (-18, 6)]
    for r, g, lab, off in zip(res, rough, labels, offsets):
        ax.annotate(lab, (r, g), xytext=off, textcoords="offset points", fontsize=7)
    ax.set(title="(c) Changing the penalty",
           xlabel=r"Residual norm $\|Am-d\|_2$", ylabel=r"Roughness $\|Lm\|_2$")
    fig.savefig(OUT / "fig08_report.png", dpi=260)
    plt.close(fig)


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
        ax.set(title=tag, xlabel="x (pitches)")
        if col == 0:
            ax.set_ylabel("y (pitches)")

    sigma = D.std(axis=0, ddof=1)
    rng = np.random.default_rng(SEED)
    peaks = []
    for _ in range(N_NOISE):
        mk = invert(A, rng.normal(0.0, sigma), LAMBDA_REF, L)
        peaks.append(np.abs(mk).max() * 1e6)
    det = float(np.percentile(peaks, 95))
    ax = fig.add_subplot(gs[0, 2])
    ax.hist(peaks, bins=26, color="#7fb3a9", edgecolor="white")
    ax.axvline(det, color="#c0392b", ls="--", lw=1.2)
    ax.annotate(f"95th percentile\n{det:.2f}", (det, ax.get_ylim()[1] * 0.72),
                xytext=(6, 0), textcoords="offset points", fontsize=7.5)
    ax.set(title="(c) Noise-only peak", xlabel=r"Peak $|m|$ ($10^{-6}$)", ylabel="Count")
    fig.savefig(OUT / "fig09_report.png", dpi=260)
    plt.close(fig)
    print(f"[report figs] detection limit check: {det:.2f} e-6 (report: 2.74)")


def main():
    A, D, meta = load()
    d = D.mean(axis=0)
    L = diff_op(A.shape[1], SIDE)
    fig08(A, d, L)
    fig09(A, D, L)
    print(f"[report figs] written to {OUT}")


if __name__ == "__main__":
    main()
