# -*- coding: utf-8 -*-
"""
L3 — two-dimensional inverse: Hall-conductivity map, resolution, detection limit,
and the pass/fail decision for the announced 3 % uniformity specification.

The forward relation is d = A m, where A is blind to a uniform offset (every row
sums to zero), so the inversion must impose a zero-mean constraint on m.

Expected files in data/:
    L3_A.npy      (128, 100) forward matrix
    L3_d.csv      repeated readings, columns like measurement_index, repeat, d_raw
    L3_meta.json  grid geometry (10 x 10), pitch / extent, reference value, spec
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "output"
OUT.mkdir(parents=True, exist_ok=True)

A_PATH = DATA / "L3_A.npy"
D_PATH = DATA / "L3_d.csv"
META_PATH = DATA / "L3_meta.json"

SPEC_DEFAULT = 3.0          # per cent
LAMBDA_SCALES = (0.01, 0.1, 1.0, 10.0, 100.0)
N_NOISE = 500
SEED = 0


def load():
    A = np.load(A_PATH)
    meta = json.loads(META_PATH.read_text(encoding="utf-8-sig"))
    raw = np.genfromtxt(D_PATH, delimiter=",", names=True)
    return A, raw, meta


def repeats_matrix(raw):
    """Return (n_readings, n_repeats) array of readings."""
    names = raw.dtype.names or ()
    idx = next((n for n in names if "index" in n.lower()), None)
    rep = next((n for n in names if "repeat" in n.lower()), None)
    val = next((n for n in names if n not in (idx, rep)), names[0])
    if idx is None or rep is None:
        return raw[val].reshape(-1, 1)
    order = np.unique(raw[idx]).astype(int)
    reps = np.unique(raw[rep]).astype(int)
    out = np.empty((len(order), len(reps)))
    for i, k in enumerate(order):
        for j, r in enumerate(reps):
            sel = (raw[idx] == k) & (raw[rep] == r)
            out[i, j] = raw[val][sel].mean()
    return out


def difference_operator(n_pix):
    """Neighbour differences on a 1-D flattened grid (no wrap-around)."""
    rows = []
    side = int(round(np.sqrt(n_pix)))
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


def solve(A, d, lam, L):
    """Regularised least squares, then project onto the zero-mean constraint.

    A @ 1 = 0, so subtracting the mean from m leaves the predictions unchanged;
    this is exactly the constraint the probe geometry demands.
    """
    n = A.shape[1]
    lhs = A.T @ A + lam * (L.T @ L) + 1e-12 * np.eye(n)
    m = np.linalg.solve(lhs, A.T @ d)
    return m - m.mean()


def reconstruct(A, d, lam):
    L = difference_operator(A.shape[1])
    return solve(A, d, lam, L)


def map_percent(m, reference):
    return (m / reference * 100.0).reshape(10, 10)


def point_spread(A, lam, side=10):
    """Impulse response -> width of the positive lobe from its second moment."""
    n = A.shape[1]
    centre = (side // 2) * side + side // 2
    e = np.zeros(n)
    e[centre] = 1.0
    m = reconstruct(A, A @ e, lam)
    m = m - m.mean()
    lobe = np.clip(m.reshape(side, side), 0, None)
    total = lobe.sum()
    yy, xx = np.mgrid[0:side, 0:side]
    if total <= 0:
        return float("nan")
    cy, cx = (lobe * yy).sum() / total, (lobe * xx).sum() / total
    var = (lobe * ((yy - cy) ** 2 + (xx - cx) ** 2)).sum() / total
    return float(2.0 * np.sqrt(var))          # 2 sigma as the resolved width


def detection_limit(A, lam, sigma, reference, rng, side=10):
    worst = []
    for _ in range(N_NOISE):
        noise = rng.normal(0.0, sigma, size=A.shape[0])
        m = reconstruct(A, noise, lam)
        worst.append(np.abs(m / reference * 100.0).max())
    return float(np.percentile(worst, 95))


def main():
    missing = [p.name for p in (A_PATH, D_PATH, META_PATH) if not p.exists()]
    if missing:
        print(f"[L3] data not found ({', '.join(missing)}) — L3 skipped, results left empty.")
        return None

    A, raw, meta = load()
    reps = repeats_matrix(raw)
    d = reps.mean(axis=1)
    sigma = float(np.mean(reps.std(axis=1, ddof=1)))
    reference = float(meta.get("reference_conductivity", meta.get("reference", d.mean())))
    spec = float(meta.get("specification_percent", SPEC_DEFAULT))
    lam0 = float(meta.get("regularisation", 1e-2))

    # regularisation sweep
    sweep = []
    for s in LAMBDA_SCALES:
        m = reconstruct(A, d, lam0 * s)
        sweep.append((s, float(np.linalg.norm(A @ m - d)), float(np.linalg.norm(np.diff(
            np.sort(m))))))
    lam = lam0

    m = reconstruct(A, d, lam)
    grid = map_percent(m, reference)
    r_max, c_max = np.unravel_index(np.argmax(grid), grid.shape)
    r_min, c_min = np.unravel_index(np.argmin(grid), grid.shape)

    pitch = float(meta.get("probe_pitch", meta.get("pitch", 0.6)))
    resolution = point_spread(A, lam)
    rng = np.random.default_rng(SEED)
    det = detection_limit(A, lam, sigma, reference, rng)

    margin = float(np.abs(grid).max() + det)
    passed = bool(margin < spec)

    result = {
        "map": grid.round(6).tolist(),
        "features": [
            {"row": int(r_max), "col": int(c_max), "deviation_pct": float(grid.max()),
             "x_pitch": float(c_max * pitch), "y_pitch": float(r_max * pitch), "sign": "+"},
            {"row": int(r_min), "col": int(c_min), "deviation_pct": float(grid.min()),
             "x_pitch": float(c_min * pitch), "y_pitch": float(r_min * pitch), "sign": "-"},
        ],
        "resolution_probe_pitch": resolution,
        "detection_limit_fraction": det / 100.0,
        "pass": passed,
        "confidence": "95% (detection limit from 500 noise-only reconstructions)",
    }
    (OUT / "L3_results.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "map"}, indent=2))
    return result


if __name__ == "__main__":
    main()
