# AI usage log — Group G11

Every substantive use of an AI tool. One row per use; keep it honest and specific.

| # | Tool | Member | Date | Purpose | What it produced | How we verified it |
|---|------|--------|------|---------|------------------|--------------------|
| 1 | ChatGPT / Claude / Copilot (name it) | Wang Boyi | 2026-11-1x | draft the trapezoidal integration loop for the band-edge count | 20 lines of NumPy | compared against Gauss–Legendre and the closed form (agreement 4e-7) |
| 2 | | Su Wenbo | | bootstrap resampling routine | loop + percentile code | checked 193/200 vs 7/200 split against an independent run with a different seed |
| 3 | | | | | | |
| 4 | | | | | | |
| 5 | | | | | | |

**Rules we applied**

- The prompt always stated the physical model and the assumptions, not just the desired code.
- Every generated block has an assertion or a limit check (units, signs, neutrality point).
- At least one member can explain every line that made it into the pipeline.
- Errors that the AI introduced are recorded in i_error_ledger.md, not here.
