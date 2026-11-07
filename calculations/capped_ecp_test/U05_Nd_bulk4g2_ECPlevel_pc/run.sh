#!/bin/bash
# ECP-level spin-free test, bare point charges: step0 free ion -> step1 in field -> step2 SA-CASSCF + AILFT (orca.out)
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
for s in step0 step1; do
  "$ORCA_BIN" $s.inp > $s.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" $s.out || { echo ${s}_FAIL; tail -25 $s.out; exit 1; }
done
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
