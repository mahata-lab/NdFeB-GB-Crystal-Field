# Single-cell embedding (comparison only)

The rare-earth ion embedded in the point charges of ONE unwrapped supercell (67 bulk / 135 GB charges, Fe charge 0), which is
not charge neutral (net +9 e bulk, +21 e boundary). This embedding reverses the sign of A20<r2> at the boundary site and gives a
spurious near-basal axis (Table 3, Fig. 3a open symbols, Supplementary Table S5). The point-charge files are kept exactly as
used; they are not produced by embedding/build_inputs.py. Parsed results: analysis/parsed/PUB_<Ion>_<bulk|GB>.json.

The Gd single-cell values (PUB_Gd_*) were obtained in the same embedding as step 2 of a separate susceptibility calculation
(SA-CASSCF(7,7), 1 octet + 48 sextets, read from converged orbitals); their inputs are not included because that job was set up
differently. They enter only the single-cell rows of Table S5 and Fig. 3a.
