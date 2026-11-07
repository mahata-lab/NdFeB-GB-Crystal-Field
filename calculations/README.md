# Calculations

Every ORCA job of the paper, one folder per job: inputs, point-charge embedding, driver script and provenance.
Output files are not included; the quantities extracted from each `orca.out` are in `../analysis/parsed/<job>.json`.
`jobs.csv` lists every job with its site, embedding, method, status and the table or figure that uses it.

| Folder | Jobs | Content | Used in |
|---|---|---|---|
| `converged_sphere/` | R01–R28 | Ce–Dy at bulk 4g (atom 2) and GB class A (atom 13), R = 18 and 25 Å; R25–R28 alternative charge partition (B −1, Fe −5/14) for Nd, Dy | Tables 2–3, Figs 2–5 |
| `matched_sites/` | T01–T12 | Ce–Dy at bulk 4f (atom 5) and GB class B (atom 9), R = 25 Å | Table 2, Figs 2–4 |
| `nevpt2_production/` | T13–T16 | Nd, Dy at bulk 4g / GB class A with SC-NEVPT2 | Table 3, Table S5 |
| `ce_active_space_tests/` | T17–T22 | Ce restart (T17) and CASCI on free-ion orbitals (T18–T22) | Table S6 |
| `capped_ecp_test/` | U01–U16 | Nd, Dy at the four matched sites; bare point charges (`_pc`) vs capped ECPs on all Fe/RE within 4.0 Å (`_cecp`) | Table 3, Table S8 |
| `boundary_chemistry/` | V01–V08 | Nd, Dy at GB classes A and B with the nearest Fe replaced by B or the nearest RE by Zr (unrelaxed) | Table 4, Table S9 |
| `midgrain_sites/` | W01–W04 | Nd, Dy at GB classes C (atom 2) and D (atom 6) | Table 4, Table S9 |
| `single_cell_embedding/` | 10 | Embedding in the charges of one non-neutral supercell (comparison) | Table 3, Fig 3a, Table S5 |
| `single_cell_basis_test/` | C07, C08 | SARC2-DKH-TZVP in the single-cell embedding | Table S5 |

Four jobs did not reach a pure 4f<sup>n</sup> active space and are kept to document this (R04, T01, T17, V01; see `NOTE.md` in
each folder). All other jobs terminated normally with seven pure 4f active orbitals.

## Files in a job folder

| File | Purpose |
|---|---|
| `sphere.pc` | Point charges (ORCA format: number of charges, then `q x y z` in Å) |
| `step0.inp` | Free ion, no field: clean atomic 4f orbitals (identical for all jobs of one ion and level) |
| `step1.inp` | Ion in the point-charge field, SA-CASSCF from the step0 orbitals (`ActConstraints 1`) |
| `step2.inp` | Final SA-CASSCF + QDPT spin–orbit coupling + AILFT + ground-doublet *g* tensor → `orca.out` |
| `step0a.inp`, `step0_rohf.inp`, `pick_f_rotations.py` | Capped-ECP arm only: ROHF guess in the same environment, rotation of the 4f orbitals into the active window, SA-CASSCF guess |
| `run.sh` | Runs the stages in order; requires `ORCA_BIN` (full path of the `orca` binary) |
| `provenance.json` | Centre, radius, charges, residual charge and its neutralisation, number of charges, nearest-charge distance |

`%pal nprocs` and `%maxcore` in the inputs are those of the original runs (16 cores, 14 GB per core for the DKH jobs; 12 cores for
the ECP jobs); adjust them to your machine. A job takes about 10–40 min on 12–16 cores (longer with SC-NEVPT2); stage 0 of Ce can need ~9 GB per process.

The four capped-ECP jobs U02, U08, U10 and U12 were run with earlier versions of the step-0 guess (default guess, or ROHF without
rotation); the scripted ROHF + rotation route in `run.sh` converges to the same seven-4f active space.
