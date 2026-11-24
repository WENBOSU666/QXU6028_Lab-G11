# L3 brief note — two-dimensional inverse

## Question

Where does the Hall conductivity of the G11 wafer deviate from the reference, what feature
size can the probe resolve, and does the wafer pass the announced 3.0 % uniformity
specification?

## Data

| File | Shape | Meaning |
|---|---|---|
| `L3_A.npy` | 128 × 100 | forward matrix, already row-centred (row sums ≈ 1e-16) |
| `L3_d.csv` | 3 × 128 | three repeats of the 128 readings, no header |
| `L3_meta.json` | — | 10 × 10 grid over a 6s × 6s region, s = probe pitch |

Because every row of `A` sums to zero, the probe cannot see a uniform offset: adding a
constant to every pixel leaves all 128 readings unchanged. The inversion therefore has to
impose a **zero-mean constraint**, otherwise the solution is dominated by an arbitrary
constant. `rank(A) = 92`, so eight directions of the 100-pixel map are unobservable.

## Method

1. Average the three repeats at each channel; noise σ per channel = std over its three repeats (mean 1.61 × 10⁻⁷, max 4.63 × 10⁻⁷).
2. Solve the regularised normal equations `(AᵀA + λ LᵀL) m = Aᵀ d` with `L` the first-order
   neighbour-difference operator on the 10 × 10 grid, then subtract the mean of `m`
   (exactly the zero-mean constraint, since `A @ 1 = 0`).
3. Reference smoothness η₀ = 10⁻² · ‖A‖₂² / ‖L‖₂² = 1.331 × 10⁻³ (methods guide); the
   decision is checked at 0.01×, 0.1×, 1×, 10× and 100× this value.
4. Resolution from a zero-mean unit-impulse test map: reconstruct with the same η and solver,
   then FWHM = 2.355 σ with σ² = Σ(w r²)/(2 Σw) over the positive lobe, r measured from the
   injected pixel (methods guide, p. 12).
5. Detection limit from 500 noise-only reconstructions (seed 20260924): the 95th percentile of
   the largest absolute deviation over the map, evaluated at the same η as the reported peak.

## Results

| Quantity | Value |
|---|---|
| maximum deviation | **+2.9094 %** at (−0.9, −2.1) probe pitches |
| minimum deviation | **−1.7557 %** at (+2.7, +2.7) probe pitches |
| spatial resolution | **1.444 probe pitches** (centre 1.444, off-centre 1.304) |
| detection limit (95 %) | **0.000274 %** (per-channel noise, 500 draws, seed 20260924) |
| noise-adjusted margin | **2.9097 %** |
| decision | **pass** (2.9097 % < 3.0 %) |

Across the five smoothness strengths the peak position never moves and the amplitude runs
from 2.9540 % (weakest) to 2.5263 % (strongest), so the pass decision is stable; the weakest
smoothing gives the conservative bound 2.9558 % after adding the detection limit.

Caveat: smoothing biases narrow features — a genuine +4.00 % single-pixel feature
reconstructs at about +1.97 % — so the reported peak is a lower bound on the true amplitude,
and the region around (−0.9, −2.1) would be rescanned at a finer pitch before fabrication.

## Figures

The report figures live in docs/figures/report/ and are produced by src/report_figures.py:

| Report figure | File | Panel |
|---|---|---|
| Figure 8(a) | ig08_report.png | zero-mean non-uniformity map, max circled, min square-marked |
| Figure 8(b) | ig08_report.png | singular values of the forward matrix (rank 92) |
| Figure 8(c) | ig08_report.png | roughness versus residual norm across the five penalty strengths |
| Figure 9(a),(b) | ig09_report.png | reconstructed unit impulse, centre and off-centre |
| Figure 9(c) | ig09_report.png | 500 noise-only reconstructions and the 95th percentile |

The five analysis panels (ig_L3_map, ig_L3_spectrum, ig_L3_sweep, ig_L3_psf,
ig_L3_noise) are written to output/ on every run of 
eproduce.py and are not tracked,
because the report figures above already contain the same information.
