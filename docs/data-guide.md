# Data dictionary and paper mapping

All paths below are relative to the repository root. `calculations/jobs.csv` maps the complete coefficient archive to its source datasets and scientific scope. It includes incomplete jobs without numerical records. `analysis/source_manifest.json` records repository paths and SHA-256 hashes of original outputs; original-source strings describe archive provenance and are not required local paths.

| Figure/table | Numerical source |
|---|---|
| Figure 1; Tables S1-S3 | Two CONTCAR structures and `structures/site_coordination.csv`; retained structural image |
| Main Tables 1-2; Figure 2 | `analysis/parsed/FINITE_*.json`; `finite_cell_descriptors.csv` |
| Figure 3 | Zero-field 2 K tensors parsed from `calculations/finite_cell_magnetic/*/step2_suscep.out.gz`; `finite_cell_susceptibility_tensors.csv` |
| Figure 4 | Zero-field 1-40 K trace averages from the same outputs; `finite_cell_chiT.csv` (480 rows) |
| Figure 5; Table S10 | Six full `analysis/orbitals/*.cube.gz`; `orbital_extent.csv`, `orbital_slice_provenance.json` |
| Figure S1; Tables S6-S7 | FINITE/R/T/N records, with model and parent sites kept separate |
| Figure S2; Table S5 | Exact rank invariants from the same AILFT matrices |
| Table S4 | FINITE many-electron QDPT-SOC spectra |
| Table S8 | U01-U16 matched-level bare-charge/capped-ECP arms, with R/T DKH reference |
| Table S9 | Direct finite-cell axis angles and the conditional Stoner-Wohlfarth expression |

## Numerical fields

- `split`: spread of seven one-electron AILFT energies, cm^-1.
- `orb_d1`: epsilon2-epsilon1; `orb_d2`: epsilon3-epsilon1. The main paper's second adjacent orbital gap Delta_orb is **orb_d2-orb_d1**, not orb_d2 alone.
- `E1`: first excited Kramers doublet for Ce/Nd/Gd/Dy, or next level above the lowest non-Kramers pair for Pr/Tb; `tun`: lowest-pair splitting. `Jspan`: assigned ground-J manifold span. All cm^-1.
- `B20`, `B40`, `B60` and `Re_Bkq`, `Im_Bkq`: Wybourne coefficients in cm^-1, cell frame z parallel c and x parallel a. Negative q follows B(k,-q)=(-1)^q conjugate(B(k,q)). `A20r2=B20/2`.
- `S2/S4/S6`: sqrt(sum_q |B(k,q)|^2/(2k+1)), cm^-1. Exact angular-rank invariants; not ligand-shell or radial decompositions.
- `g_angle`: direct ground-doublet largest-g axis tilt from c when available. It is not the main Figure 3 bulk-to-boundary susceptibility-axis angle. The main latter angle is `bulk_GB_axis_angle_deg`.
- `chiT`: molar orientational average T Tr(chi)/3 in cm3 K/mol, per rare-earth ion. Tensor entries/eigenvalues in the 2 K table are likewise chi*T in cm3 K/mol. Unit factor T does not change their eigenvectors. Tensor axis signs are arbitrary.
- `grid_norm`: integral of |psi|^2 over the full cube before normalization. `orbital_r_rms_A`: full 3D one-electron eigenorbital rms radius in Å. Projected density has units Å^-2 and uses one common 0–7 Å^-2 range with a square-root display coordinate for all panels; colorbar ticks remain in physical units.
- `conditional_Hsw_over_HA`: dimensionless Stoner-Wohlfarth loss-of-metastability field; not a measured/predicted coercivity or spatial nucleation field.

`cf_all_runs.csv` and `wybourne_Bkq_all_runs.csv` contain separate rows for CASSCF and NEVPT2 where available. The complete coefficient CSV retains all 124 run/block rows of the current supplementary data, including historical records. It is not embedded in the supplementary PDF. Do not pool different models, method blocks or duplicate Gd aliases as independent observations.

`finite_cell_run_reconciliation.csv` compares FINITE with the separate PUB outputs. For example, Ce bulk orbital spread is 2589.4 cm^-1 in the magnetic-property rerun and 2594.6 cm^-1 in the separate production output. These are different saved calculations, not a rounding correction.

The Figure 4 data are one-third of the trace of the explicitly zero-field printed chiT tensors, retaining the original 480 curve values to their saved rounding. `zero_field_response.py` selects the field headers explicitly and does not read the separate scalar susceptibility column.

The [axis-source reconciliation](../analysis/tables/axis_source_reconciliation.csv) records the reproduced submitted bars, the corrected original partial-g direction indicator and the direct susceptibility-axis result. The auxiliary-reference overwrite is explained in [model scope](model-scope.md). The Gd values are diagnostics only, not an assigned orbital easy axis.
