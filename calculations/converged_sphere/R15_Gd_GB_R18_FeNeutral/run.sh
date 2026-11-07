#!/bin/bash
# step0 free ion -> step1 ion in the point-charge field -> step2 SA-CASSCF + SOC + AILFT (result: orca.out)
# step0 is identical for every job of a given ion; if step0.gbw is supplied it is reused.
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
if [ ! -f step0.gbw ]; then
  "$ORCA_BIN" step0.inp > step0.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" step0.out || { echo STEP0_FAIL; tail -25 step0.out; exit 1; }
fi
echo STEP0_OK
"$ORCA_BIN" step1.inp > step1.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step1.out || { echo STEP1_FAIL; tail -25 step1.out; exit 1; }
echo STEP1_OK
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
grep -q "active space contains 7 f Orbitals : OK" orca.out && echo "STEP2_OK (7 pure 4f active orbitals)" || echo "STEP2_WARN: check active space"
