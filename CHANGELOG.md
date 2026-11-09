# Changes

## Figure 5 label and contour update

Added ion/environment/MO27 panel titles, labeled 0.5 Å spatial ticks, 0.25 Å minor ticks and thin black isodensity contours at shared absolute levels 0.25, 0.75 and 2.0 Å^-2. Updated the related caption and explanations only. Source numerical tables and other figure topics are unchanged.

## Figure 5 contrast update

Retained the current cube data, projection, six panels and radius annotations. Improved weak-density visibility with a shared square-root color mapping over 0–7 Å^-2, with physical-density colorbar ticks. Updated the caption/discussion, response, supplement, provenance and reproduction script; no numerical tables or simulations changed.

## 2026-11-07 — current revision data

- Restored the finite-cell main comparison and all five original figure topics, including the six-panel susceptibility plot.
- Added the twelve saved magnetic-property inputs and complete losslessly compressed outputs; distinguished FINITE reruns from separate PUB production datasets.
- Recomputed magnetic response directions from saved zero-field 2 K susceptibility tensors, omitting a Gd orbital easy-axis assignment.
- Replaced separately normalized orbital slices with six full-cube density projections on one absolute scale; supplied the six compressed cubes and full 3D extent metrics.
- Added the focused supplement: Figures S1-S2, Tables S1-S10 and the complete 124-row complex Wybourne CSV.
- Indexed matched parent sites, embedding/radius/charge-partition/NEVPT2/cECP controls, active-space failures and historical explorations without treating their models as interchangeable.
- Updated the title and manuscript status; removed unsupported statements that all original finite-cell results were spurious, that neutrality uniquely validates the model, or that orbital gaps directly predict coercivity.
- Preserved older figures, tables and scripts under `archive/`; retained calculation identifiers and input paths.
- Added a portable reproduction entry point and source/value verification. No new ORCA simulations were run, and this local update was not committed or pushed.
