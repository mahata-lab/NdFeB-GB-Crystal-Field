#!/usr/bin/env bash
set -euo pipefail
: "${ORCA_BIN:?Set ORCA_BIN to the full path of the ORCA executable}"
# This archived magnetic-property input reads the saved production orbitals.
test -s step2.gbw || { echo "Supply the converged finite-cell step2.gbw; see ../README.md" >&2; exit 1; }
"$ORCA_BIN" step2_suscep.inp > step2_suscep.out
