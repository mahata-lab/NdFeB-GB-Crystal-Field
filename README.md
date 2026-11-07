# Grain-boundary 4f crystal fields in Nd₂Fe₁₄B — inputs, embeddings and analysis

Data and code for

> A. Mahata and M. Zakotnik, *Grain Boundary Distortion Strengthens the Local 4f Crystal Field and Separates Heavy from
> Light Lanthanides in Nd₂Fe₁₄B*, Materialia (2026).

The repository contains everything needed to rebuild the point-charge embeddings, rerun every ORCA calculation, and
regenerate all tables and figures of the paper from the extracted results.

## Contents

```
structures/     relaxed DFT supercells (VASP CONTCAR) and the Fe coordination of every rare-earth site
                  bulk_Nd8Fe56B4_CONTCAR          bulk Nd2Fe14B, 68 atoms
                  GB90twist_Nd16Fe112B8_CONTCAR   90° [001] twist grain boundary, 136 atoms
                  site_coordination.csv           (Supplementary Table S3)
embedding/      build_inputs.py      regenerates every job folder in calculations/ from structures/ (and checks them)
                pick_f_rotations.py  4f orbital selection for the capped-ECP jobs
                run_jobs.sh          runs job folders one after another on a workstation
                submit_template.slurm
calculations/   one folder per ORCA job: step*.inp, sphere.pc, run.sh, provenance.json   (see calculations/README.md)
                jobs.csv             index of all 90 jobs: site, embedding, method, status, where used in the paper
analysis/       parsed/   results extracted from every orca.out (AILFT matrix, Slater-Condon parameters, zeta,
                          spin-orbit levels, Kramers-doublet <L>,<S>), one JSON per job
                scripts/  parsing, crystal-field analysis, tables, figures
                tables/   all derived tables (CSV)
figures/        Figures 1-5 of the paper
```

## Requirements

* ORCA 6.1.1 (https://orcaforum.kofo.mpg.de) with the matching OpenMPI, for the quantum-chemical calculations.
* Python ≥ 3.9 with numpy, scipy (≥ 1.15), pandas and matplotlib (`pip install -r requirements.txt`) for everything else.

## Reproducing the paper

**1. Embeddings and inputs.** `python3 embedding/build_inputs.py --check` rebuilds all 78 generated job folders from the two
CONTCAR files and compares them with `calculations/` (point charges to 10⁻⁶ e/Å, inputs line by line). To build fresh
folders: `python3 embedding/build_inputs.py OUTDIR [job prefixes]`.

**2. ORCA calculations.** In a job folder, `ORCA_BIN=/full/path/to/orca bash run.sh`, or for several folders
`bash embedding/run_jobs.sh calculations/matched_sites/T0*`. Each job runs three stages (free ion → ion in the field →
SA-CASSCF + QDPT spin–orbit coupling + AILFT + g tensor) and writes `orca.out`.

**3. Extraction.** `python3 analysis/scripts/orca_parse.py <job>/orca.out <job> <ion> analysis/parsed/<job>.json`.
The JSON files of all published calculations are included, so steps 4–5 can be run directly.

**4. Tables** (from `analysis/`):

```
python3 scripts/cf_tables.py parsed tables/cf_all_runs.csv            # all descriptors, one row per calculation
python3 scripts/derived_tables.py tables/cf_all_runs.csv tables/      # Table 2 (matched pairs), Table 3 / S5 (sensitivity)
python3 scripts/uvw_tables.py tables/cf_all_runs.csv tables/uvw_tests.csv   # capped-ECP, substitution, mid-grain (Tables 3, 4)
python3 scripts/export_bkq.py parsed tables/wybourne_Bkq_all_runs.csv # complete Wybourne parameter sets (Table S4)
python3 scripts/easy_axis_lfci.py parsed tables/easy_axis_lfci.csv $(cut -d, -f1 tables/easy_axis_lfci.csv | tail -n +2)
python3 scripts/nucleation_1d.py tables/nucleation_1d.csv              # Table S7
python3 scripts/pointcharge_sites.py ../structures/bulk_Nd8Fe56B4_CONTCAR ../structures/GB90twist_Nd16Fe112B8_CONTCAR tables 25
python3 scripts/make_supplementary.py . ../structures/bulk_Nd8Fe56B4_CONTCAR ../structures/GB90twist_Nd16Fe112B8_CONTCAR ../structures/site_coordination.csv si_tables.tex
```

**5. Figures:** `python3 scripts/make_figures.py . ../figures` (Figures 2–5; Figure 1 is the structure rendering).

Running steps 4–5 on the included JSON files reproduces the tables in `analysis/tables/`, the Supplementary tables and
Figures 2–5 exactly.

## Conventions

* **Embedding.** Point charges RE +3, B −3, Fe −3/14 (formula-unit neutral; alternative partition B −1, Fe −5/14). All periodic
  images within R of the ion are included; the residual charge of the truncated sphere is spread over the charges in the outer
  2 Å, so every embedding is exactly neutral. Production radius R = 25 Å. The ion is a bare quantum ion (no explicit ligands).
* **Sites.** One-based atom indices of the CONTCAR files. Bulk atoms 1–4 are Nd 4g (x,−x,0), atoms 5–8 Nd 4f (x,x,0)
  (Herbst et al. coordinates). Grain-boundary classes: A {5,8,13,16} and B {1,4,9,12} on the boundary planes (images of 4f and 4g),
  C {2,3,10,11} and D {6,7,14,15} midway between the planes. Matched pairs: GB atom 13 (A) ↔ bulk atom 5 (4f);
  GB atom 9 (B) ↔ bulk atom 2 (4g).
* **Descriptors** (cm⁻¹). Total 4f splitting = spread of the seven AILFT orbital energies; Wybourne parameters B<sup>k</sup><sub>q</sub>
  from an exact decomposition of the AILFT matrix (cell frame, z ∥ c); A₂⁰⟨r²⟩ = B²₀/2; S<sub>k</sub> = [Σ<sub>q</sub>|B<sup>k</sup><sub>q</sub>|²/(2k+1)]<sup>1/2</sup>;
  E₁ = first excitation of the ground J multiplet from the ab initio spin–orbit spectrum. ORCA prints the AILFT matrix in
  Hartree; `orca_parse.py` converts to cm⁻¹.

## Licence

Code: MIT (see `LICENSE`). Data (structures, inputs, parsed results, tables, figures): CC BY 4.0. Please cite the paper when
using these data.
