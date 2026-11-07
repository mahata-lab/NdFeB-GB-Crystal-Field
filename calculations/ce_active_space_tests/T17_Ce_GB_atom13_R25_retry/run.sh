#!/bin/bash
# step1 starts from supplied orbitals (see NOTE.md), then step2 (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
"$ORCA_BIN" step1.inp > step1.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step1.out || { echo STEP1_FAIL; tail -25 step1.out; exit 1; }
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
