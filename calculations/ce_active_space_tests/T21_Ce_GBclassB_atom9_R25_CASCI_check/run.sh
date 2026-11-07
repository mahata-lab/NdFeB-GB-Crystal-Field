#!/bin/bash
# Ce CASCI variant: step0 free ion, then step2 = SA-CASCI (NoIter) on the free-ion 4f orbitals + SOC + AILFT (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
if [ ! -f step0.gbw ]; then
  "$ORCA_BIN" step0.inp > step0.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" step0.out || { echo STEP0_FAIL; tail -25 step0.out; exit 1; }
fi
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
