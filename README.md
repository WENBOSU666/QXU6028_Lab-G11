# QXU6028 — Group G11

**Carrier counting and conductivity uniformity in graphene on sapphire**

Reproducible analysis for the QXU6028 virtual lab: L1 forward carrier counting, L2
zero-dimensional inverse (DoS class, band gap, Fermi-level trajectory), L3 two-dimensional
inverse (uniformity map, resolution, detection limit, pass/fail).

## Quick start

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
python reproduce.py
```

The command reads the raw files in `data/`, runs the three stages in order and writes
`output/results.json` and `results.json` (schema `QXU6028-results-v1`). Runtime is about
10 s for L1 and L2.

## Repository layout

```
.
├── reproduce.py            one command: L1 -> L2 -> L3 -> merged results.json
├── requirements.txt
├── results.json            final numbers in the submitted schema (regenerated)
├── src/
│   ├── l1_forward.py       band-edge vs Fermi-referenced counting, correction factor
│   ├── l2_inverse.py       DoS class, band gap, Fermi-level trajectory, bootstrap
│   └── l3_map.py           zero-mean regularised inversion, PSF, detection limit
├── data/                   supplied datasets for Group G11
├── output/                 figures and per-stage JSON (git-ignored)
└── docs/
    ├── figures/            figure files used in the report
    ├── L1_brief_note.md, L2_brief_note.docx
    ├── ai_usage_log.md
    ├── ai_error_ledger.md
    └── contribution_statement.docx / .pdf
```

## Results

| Stage | Result | Where |
|---|---|---|
| L1 | correction factor **C = 3.953** at E_F = 0.15550 eV; n_F = 4.6127 × 10¹¹ cm⁻²; mobility **47,337 cm² V⁻¹ s⁻¹** | `output/fig_L1.png` |
| L2 | **hard_gap** model; **E_g = 3.25 k_BT (84.0 meV)**; E_F from −77.9 to +135.7 meV; RMS 0.0987 → 0.0126 | `output/fig_L2_*.png` |
| L3 | map max **+2.9094 %** at (-0.9, -2.1), min -1.7557 %; resolution **1.74 probe pitches**; detection limit 0.00022 %; **pass** (margin 2.9096 % < 3.0 %) | output/fig_L3_*.png |

## Reproducibility notes

- Both stages read only from `data/` and write only to `output/`; there are no absolute paths.
- Random procedures are seeded (`reproduce.py --seed 0`): the L2 bootstrap uses a fixed
  generator seed and the L3 noise study uses `default_rng(0)`.
- `results.json` records `pipeline_commit`, the commit the numbers were produced from.
- A clean-clone test is documented in the report: clone into an empty folder and run
  `python reproduce.py`.

## L3 status

L3 is complete. src/l3_map.py reconstructs the zero-mean map, locates the extrema,
estimates the resolution from the point-spread function and the detection limit from 500
noise-only reconstructions, then decides pass/fail against the 3.0 % specification.
See docs/L3_brief_note.md for the method and output/fig_L3_*.png for the figures.

## Team

| Member | ID | Phases |
|---|---|---|
| Wang Boyi | 2023304384 | L1 + L2 + L3 |
| Su Wenbo | 2023304316 | L1 + L2 + L3 |
| Lin Zhiyuan | 2023304268 | L1 + L2 |
| Zhai Xiangru | 2023304281 | L2 + L3 |

Signed record: `docs/contribution_statement.pdf`.
