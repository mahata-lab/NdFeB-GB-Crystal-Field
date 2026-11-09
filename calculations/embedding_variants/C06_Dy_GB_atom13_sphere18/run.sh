#!/bin/bash
# three-step driver: atomic guess (step0) -> field (step1) -> full+AILFT (step2)
set -e
echo "=== STEP0 : free-ion atomic 4f ==="
"$ORCA_BIN" step0.inp > step0.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step0.out || { echo STEP0_FAIL; tail -25 step0.out; exit 1; }
echo STEP0_OK

echo "=== STEP1 : field, reads atomic guess ==="
"$ORCA_BIN" step1.inp > step1.out 2>&1 || true
grep -q "ORCA TERMINATED NORMALLY" step1.out || { echo STEP1_FAIL; tail -25 step1.out; exit 1; }
grep -qi "contains 0 f Orbitals\|f Orbitals : not OK" step1.out && echo "STEP1_WARN: impure f (continuing)"
echo STEP1_OK

echo "=== STEP2 : full SA-CASSCF + SOC + AILFT ==="
"$ORCA_BIN" step2.inp > orca.out 2>&1 || true
if grep -q "ORCA TERMINATED NORMALLY" orca.out; then
  echo STEP2_OK
  echo "==== AILFT / f-orbital check ===="
  grep -i "f Orbitals\|AILFT\|SOC\|Stevens\|not OK" orca.out | head -20
else
  echo STEP2_FAIL
  grep -i "f Orbitals\|not OK\|AILFT" orca.out | tail -10
  tail -20 orca.out
  exit 1
fi
