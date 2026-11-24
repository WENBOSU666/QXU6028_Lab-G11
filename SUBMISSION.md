# Submission map — QXU6028 Group G11

Everything the brief asks for, and where to find it.

| # | Required item | Where |
|---|---|---|
| 1 | Group report, ≤ 2,500 words + figures, with a ≤ 200-word Transfer section | submitted as PDF through QMplus (source: `docs/` copy kept with the group) |
| 2 | Code repository, open to the markers, full commit history, one-command rerun | this repository — `python reproduce.py` |
| 3 | `results.json` in the supplied schema | `results.json` at the repository root (also uploaded separately) |
| 4 | AI error ledger | `docs/ai_error_ledger.md` |
| 5 | AI usage log | `docs/ai_usage_log.md` |
| 6 | Signed contribution statement | `docs/contribution_statement.pdf` (editable: `.docx`) |
| 7 | 150–200-word critical self-assessment | submitted with the report |

## One command

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
python reproduce.py
```

Runtime ≈ 20 s. It writes `output/results.json` and `results.json`, eleven analysis figures
and the nine report figures (`output/report/fig01_report.png` … `fig09_report.png`).
No absolute paths, no manual steps, all random steps seeded (`--seed 0`).

## What is generated vs. supplied

| Supplied (never modified) | Generated (safe to delete, regenerated) |
|---|---|
| `data/L1_*.csv/json`, `L2_*.csv/json`, `L3_*.npy/csv/json` | `output/**`, `results.json` |

The datasets were issued to Group G11 only and have not been shared with any other group.

## Result summary

| Stage | Result |
|---|---|
| L1 | correction factor **C = 3.953**; n_F = 4.6127 × 10¹¹ cm⁻²; mobility **47,337 cm² V⁻¹ s⁻¹** |
| L2 | **hard_gap**; E_g = **3.25 k_BT (84.0 meV)**; E_F −77.9 → +135.7 meV over 51 measurements |
| L3 | map max **+2.9094 %** at (−0.9, −2.1); min −1.7557 %; resolution **1.74 probe pitches**; detection limit 0.000275 %; **pass** (margin 2.9097 % < 3.0 %) |
| Blind test | the same pipeline runs unchanged on the new data set issued in session 3 |

## Repository layout

```
reproduce.py          one command: L1 → L2 → L3 → report figures → results.json
requirements.txt
results.json          final numbers, schema QXU6028-results-v1
src/
  l1_forward.py       band-edge vs Fermi-referenced counting, correction factor
  l2_inverse.py       DoS class, band gap, Fermi-level trajectory, bootstrap
  l3_map.py           zero-mean regularised inversion, PSF, detection limit, decision
  report_figures.py   the nine figures in the report's panel layout
data/                 supplied datasets (Group G11 only)
output/               generated, not tracked
docs/
  ai_usage_log.md, ai_error_ledger.md
  contribution_statement.pdf / .docx
  L1/L2/L3_brief_note, figures/report/fig01…fig09
  manifest.json, results_template.json, README_FIRST.txt (as supplied)
```

## Markers: shortest path

1. clone the repository;
2. `python reproduce.py`;
3. compare `results.json` with the submitted file — every value is reproduced;
4. open `docs/figures/report/` — those nine images are the report figures, same layout.
