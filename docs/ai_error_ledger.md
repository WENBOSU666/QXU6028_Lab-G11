# AI error ledger — Group G11

Specific mistakes produced by AI tools, how we detected them, and how we corrected them.
Two genuine physics errors are expected; generic "the tool had a bug" entries earn little credit.

## Error 1 — wrong lower integration limit in the band-edge count

- **Where**: L1 forward model (src/l1_forward.py).
- **What the AI produced**: the band-edge count integrated from E_F instead of the Dirac point,
  so the two counts became identical and C collapsed to 1.
- **How we detected it**: the L1 gate requires C ≈ 2 at E_F = 0.10 eV; our sweep gave C < 1
  at positive doping. The sign of the discrepancy pointed at the integration limits.
- **Correction**: lower limit of the band-edge integral set to the Dirac point, the Fermi-referenced
  count kept at E_F; the sweep then reproduced the reference points 1.00 / 2.599 / 5.104.
- **Evidence**: commit <hash>; Figure 1(b) before/after.

## Error 2 — normalisation applied before locating charge neutrality

- **Where**: L2 inverse (src/l2_inverse.py).
- **What the AI produced**: it normalised R_H by its global maximum before finding the zero-crossing,
  which shifted the neutrality point and broke the wings relation |R̃_H| ∝ R̃_xx².
- **How we detected it**: the L2 gate check — after normalisation R_H no longer crossed zero
  between indices 24 and 25.
- **Correction**: locate the zero-crossing first, then normalise each channel by its own maximum.
- **Evidence**: commit <hash>; Figure 5.

## Error 3 (optional) — zero-mean constraint missing in L3

- **Where**: L3 inversion (src/l3_map.py).
- **What the AI produced**: a least-squares solve of d = A m without the zero-mean constraint;
  the map came out dominated by a flat offset because every row of A sums to zero.
- **How we detected it**: the reconstruction had a large uniform component and the residual did not improve.
- **Correction**: solved the constrained problem with the mean of m fixed to zero and compared
  five regularisation strengths.
- **Evidence**: commit <hash>; Figure 8(c).
