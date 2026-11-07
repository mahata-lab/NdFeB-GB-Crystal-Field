#!/bin/bash
# ECP-level spin-free test, capped ECPs on the Fe and rare-earth neighbours within 4.0 A (orca.out).
# step0a: ROHF on the same atom list (ion + uncharged caps) with Loewdin MO populations;
# pick_f_rotations.py rotates the seven lowest MOs of mainly 4f character into the CAS window (if needed);
# step0: SA-CASSCF from those orbitals; step1/step2 in the field as for the point-charge arm.
set -e
: "${ORCA_BIN:?set ORCA_BIN to the full path of the orca binary}"
HERE="$(cd "$(dirname "$0")" && pwd)"; PICK="$HERE/pick_f_rotations.py"; [ -f "$PICK" ] || PICK="$HERE/../../../embedding/pick_f_rotations.py"
"$ORCA_BIN" step0a.inp > step0a.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step0a.out || { echo step0a_FAIL; tail -20 step0a.out; exit 1; }
python3 "$PICK" step0a.out 9 7 rotate.txt
if [ -s rotate.txt ]; then
  awk 'NR==FNR{r=r $0 "\n"; next} /^%maxcore/{print; printf "%s", r; next} {print}' rotate.txt step0_rohf.inp > step0.inp
else
  cp step0_rohf.inp step0.inp
fi
for s in step0 step1; do
  "$ORCA_BIN" $s.inp > $s.out 2>&1 || true
  grep -q "ORCA TERMINATED NORMALLY" $s.out || { echo ${s}_FAIL; tail -25 $s.out; exit 1; }
done
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" orca.out && echo STEP2_OK || { echo STEP2_FAIL; tail -25 orca.out; exit 1; }
