# Calculation index

`jobs.csv` indexes dataset IDs, input directories, extracted results, source hashes, status, scientific scope and current use. Four unsuccessful 4f-active-space jobs remain documented (R04, T01, T17, V01); they have no parsed numerical result. A normal program exit and a suitable 4f active space are different checks.

| Folder | Role in the current revision |
|---|---|
| `finite_cell_magnetic/` | Twelve main magnetic-property reruns; saved inputs, charge arrays and complete `.out.gz` outputs |
| `single_cell_embedding/` | Ten separate original production input sets, indexed as PUB datasets; not the main magnetic-property reruns |
| `finite_cell_nevpt2/` | N01-N04 finite-cell NEVPT2 checks, SI Table S6 |
| `converged_sphere/` | Historical folder name for R01-R28 neutral-sphere radius/charge checks; selected SI S6-S7 / Figs S1-S2 |
| `matched_sites/` | T01-T12 matched parent-site checks; SI S7 / Figs S1-S2 |
| `nevpt2_production/` | T13-T16 sphere-model NEVPT2 checks; SI S6; historical path name retained |
| `ce_active_space_tests/` | T17-T22 Ce restart/CASCI limitations; no substitute Ce class A main result |
| `capped_ecp_test/` | U01-U16 matched-level spin-free core-repulsion control; SI S8 |
| `single_cell_basis_test/`, `embedding_variants/` | Archived C-series basis/embedding checks |
| `boundary_chemistry/`, `midgrain_sites/`, `exploratory_sites/` | Historical V/W/S site and unrelaxed chemistry checks; not used in current figures or focused supplement |

Non-main raw outputs are not bundled, as in the earlier repository; parsed records, inputs and original-output hashes are supplied. All 12 main magnetic-property outputs are now bundled losslessly compressed. The six Figure 5 cubes are under `analysis/orbitals/`. Binary GBW restarts are not included. Historical PUB_Gd and FINITE_Gd labels refer to the same saved Gd outputs.

Input resource settings reflect the saved runs and may exceed a workstation's capacity. Archived `run.sh` scripts require Linux/WSL, ORCA and matching MPI. See `docs/reproducing.md` and `finite_cell_magnetic/README.md` before rerunning. No new ORCA calculations were performed for this repository update.
