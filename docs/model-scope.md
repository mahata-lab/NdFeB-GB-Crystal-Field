# Model definitions and limits

## Finite-cell main comparison

One quantum Ln(III) ion, no explicit ligand basis functions. Point charges represent the remaining atoms of one unwrapped cell: RE +3, Fe 0, B -3; 67 bulk and 135 interface entries, with sums +9 e and +21 e. These sums describe the embedding; they alone neither validate nor invalidate it. The quantum ion is not a point-charge entry. These arrays are not constructed as a neutral periodic Madelung field.

Bulk atom 2 is 4g; interface atom 13 is class A (4f parent). The finite-cell orbital-spread changes are approximately -36%, -13%, +3%, +23%, +24%, +26% for Ce, Pr, Nd, Gd, Tb, Dy. These are selected-site model results, not a universal interface effect. Direct 2 K susceptibility-axis rotations are about 87° for the five anisotropic ions; no Gd orbital easy axis is assigned. The original plotting reader let auxiliary `_casscf` records overwrite `_bulk`. Strict bulk/GB matching corrects that source mismatch. The archived g(OZ/SOC)-contribution direction indicator and the explicitly defined 2 K susceptibility axes agree within 0.04 degrees for the five anisotropic ions; they are distinct definitions. The separate `PUB_*` and magnetic-property `FINITE_*` outputs differ slightly; see the reconciliation table.

## Sphere sensitivity model

Periodic images inside 18 or 25 Å carry RE +3, Fe -3/14, B -3, with residual charge distributed over the outermost 2 Å. Alternative partition: Fe -5/14, B -1. The finite-cell-to-sphere comparison changes charges and truncation as well as size; only 18-to-25 Å within one prescription isolates radius. The folder name `converged_sphere` is a historical path, not a claim of unique physical convergence or a universal uncertainty bound.

Matched pairs: class A atom 13 versus bulk 4f atom 5; class B atom 9 versus bulk 4g atom 2. Ce class A is unavailable because the active orbitals acquire 6s character. CASCI checks do not justify filling that missing SA-CASSCF result. This instability is not a calculated Ce valence transition. No fitted LF-CI reconstruction supplies the current figures.

The matched neutral-sphere calculations do not retain the finite-cell light/heavy sign crossover or the large rotations. They expose sensitivity; neutrality alone does not prove spectroscopic accuracy. Capped ECPs test core repulsion at a matched def2-TZVP/central-ECP level and do not test explicit Fe-4f covalency or conduction-electron screening.

## Structure and magnetic interpretation

Atom indices are one-based in the supplied structures. Bulk 1-4 are 4g and 5-8 are 4f. Interface classes A={5,8,13,16}, B={1,4,9,12}, C={2,3,10,11}, D={6,7,14,15}. The upper half follows (-y,x,z+1/2), without an extra in-plane shift. This coherent registry model is not a representative general high-angle disorientation. Historical S08/S09 job names incorrectly say `4g_atom5`; atom 5 is 4f. Their original identifiers remain for provenance and their index uses the correct site.

Main results are SA-CASSCF. NEVPT2 is a separately labeled check. Equal total weights are assigned to multiplicity blocks; roots within each block are equally weighted. The included root spaces are not all spin-free terms for every heavy ion. The saved main magnetic inputs use the default 75 macroiterations; separate production inputs allow 200.

Orbital spread and orbital pair gaps describe the one-electron AILFT operator, not magnetic doublets. SOC excitations are many-electron quantities. Ground-doublet directions for sphere-model Ce/Nd/Dy come directly from saved magnetic matrix elements; no Pr/Tb/Gd direction is inferred for those sphere-model tables.

The main susceptibility curves use the original zero-static-field tensor-trace response, averaged as T Tr(chi)/3 per mole of rare-earth ions; the separate scalar d2E/dB2 column is not substituted. No Fe exchange is included. Figure 5 uses full orbital-27 cube projections with a shared 0–7 Å^-2 range and square-root display mapping (gamma=1/2), explicitly labeled with physical-density ticks. Individual-eigenorbital radii cannot establish the total correlated-shell radius or a unique radial mechanism. Table S9 is a dimensionless, conditional isolated uniaxial Stoner-Wohlfarth loss-of-metastability estimate; no K is inferred from orbital gaps and no macroscopic coercivity is predicted.
