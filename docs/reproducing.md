# Reproduction

## Current analysis, without ORCA

From the repository root, with Python >=3.10 and the dependencies in `requirements.txt`:

```bash
python analysis/scripts/reproduce.py --output build
python analysis/scripts/verify_repository.py
```

`build/analysis/` receives descriptors, coefficients, extracted finite-cell response, run reconciliation, conditional switching estimates and cube integrals. `build/figures/` receives Figures 1-5 and S1-S2. `build/supplementary/` receives the current supplementary body and regenerated numerical tables S4-S10; structural Tables S1-S3 use the retained full-coordinate template. The two structure files are checked against their source hashes by verification.

Figure 1 is copied from the supplied retained rendering. Figure 5 is regenerated from all voxels in each lossless cube; it is not a reconstruction from the PNG. The plotting code uses Cambria where installed and bundled DejaVu Serif otherwise; fonts can affect figure bytes and text width. Numerical checks compare values rather than PNG hashes. PNG previews are rendered from the generated vector PDFs.

## Individual operations

```bash
python analysis/scripts/cf_tables.py analysis/parsed build/cf_all_runs.csv
python analysis/scripts/export_bkq.py analysis/parsed build/wybourne_Bkq_all_runs.csv
python analysis/scripts/orca_parse.py calculations/finite_cell_magnetic/FINITE_Ce_bulk/step2_suscep.out.gz FINITE_Ce_bulk Ce build/Ce.json
```

Create the output directory before using individual operations. The parser accepts plain `.out` and lossless `.out.gz`. Full magnetic-property outputs are supplied for all twelve FINITE datasets. Other datasets are supplied as extracted JSON plus original-output hashes and input files; their raw outputs are not bundled. PUB_Gd aliases use the same output as FINITE_Gd.

## Input regeneration

```bash
python embedding/build_inputs.py --check
python embedding/build_inputs.py build/generated-inputs R10 U04
```

This checks the 78 sphere/Ce/capped-ECP/chemistry/mid-grain folders supported by the existing generator. It does not generate the original finite-cell magnetic-property or production inputs, N/C/S archived inputs, or binary restarts. The check ignores comments and resource settings when comparing inputs, and compares point charges to 2e-6 in the printed e/Å fields. Historical folder names are retained to keep the generator and job paths stable.

## Optional ORCA reruns

The archived jobs use ORCA 6.1.1 and matching MPI. Scripts run on Linux/WSL. Set `ORCA_BIN` to the full path of the ORCA executable, then run a chosen job's `run.sh`; adjust `%pal`/`%maxcore` to available resources. `embedding/run_jobs.sh` supports sequential scratch runs. No ORCA invocation is part of the reproduction command above.

Most generated jobs prepare free-ion orbitals, optimize them in the field and then perform the final SA-CASSCF/SOC/AILFT step. The **magnetic-property reruns** instead read a saved `step2.gbw`; these GBW binaries are not supplied. See `calculations/finite_cell_magnetic/README.md`. Capped-ECP jobs use a distinct ROHF/rotation guess route, as documented in their scripts. Failed Ce and substitution active spaces remain documented; do not treat a run exit code as validation of 4f character.

## Supplementary document

The supplied `supplementary/4f-supplementary.pdf` is the current checked paper asset. To compile the source, run pdfLaTeX twice from `supplementary/`, keeping `si_tables.tex` and the sibling `figures/` directory available. The generated `build/supplementary/` has the same relative layout. `wybourne_Bkq_all_runs.csv` is supplied separately in `analysis/tables/`; it is not a TeX input or part of the PDF. Manuscript drafts and reviewer correspondence are not included in this public data/code package.
