#!/bin/bash
# run_jobs.sh -- run ORCA job folders one after another on a Linux/WSL workstation.
#
#   export ORCA_BIN=/full/path/to/orca            # ORCA 6.1 with matching OpenMPI
#   bash embedding/run_jobs.sh calculations/matched_sites/T03_* calculations/converged_sphere/R10_*
#   NPROCS=8 MAXCORE=3000 bash embedding/run_jobs.sh calculations/midgrain_sites/*
#
# Each job runs in a scratch copy ($SCRATCH or ~/orca_scratch/<job>); %pal and %maxcore are set to NPROCS/MAXCORE
# (default: physical cores up to 12, 70% of free memory). orca.out and the stage outputs are copied back into the
# job folder; finished jobs (orca.out with "ORCA TERMINATED NORMALLY") are skipped, so the script can be restarted.
set -u
SCR="${SCRATCH:-$HOME/orca_scratch}"
ORCA_BIN="${ORCA_BIN:-$(command -v orca)}"
[ -x "$ORCA_BIN" ] || { echo "ERROR: set ORCA_BIN to the full path of the orca binary"; exit 1; }
ORCA_BIN="$(readlink -f "$ORCA_BIN")"
export PATH="$(dirname "$ORCA_BIN"):$PATH"
export LD_LIBRARY_PATH="$(dirname "$ORCA_BIN"):${LD_LIBRARY_PATH:-}"
NCORE=$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l); [ "$NCORE" -ge 1 ] 2>/dev/null || NCORE=$(( $(nproc) / 2 ))
DEF=$NCORE; [ "$DEF" -gt 12 ] && DEF=12
NPROCS=${NPROCS:-$DEF}
command -v mpirun >/dev/null 2>&1 || { echo "WARNING: mpirun not found -> serial run (NPROCS=1)"; NPROCS=1; }
MEM_MB=$(awk '/MemAvailable/ {print int($2/1024)}' /proc/meminfo)
MAXCORE=${MAXCORE:-$(( MEM_MB * 70 / 100 / NPROCS ))}; [ "$MAXCORE" -gt 14000 ] && MAXCORE=14000
export OMP_NUM_THREADS=1 OMPI_MCA_rmaps_base_oversubscribe=1
[ $# -eq 0 ] && { sed -n '2,9p' "$0"; exit 1; }
echo "ORCA $ORCA_BIN | NPROCS $NPROCS | MAXCORE $MAXCORE MB | scratch $SCR"

for J in "$@"; do
  J="$(readlink -f "$J")"; job="$(basename "$J")"
  [ -f "$J/run.sh" ] || { echo "[skip] $job (no run.sh)"; continue; }
  if [ -f "$J/orca.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$J/orca.out"; then echo "[done] $job"; continue; fi
  W="$SCR/$job"; rm -rf "$W"; mkdir -p "$W"
  cp "$J"/*.inp "$J"/*.pc "$J"/run.sh "$W"/; cp "$J"/*.py "$J"/*.gbw "$W"/ 2>/dev/null || true
  sed -i "s/%pal nprocs [0-9]* end/%pal nprocs $NPROCS end/; s/%maxcore [0-9]*/%maxcore $MAXCORE/" "$W"/*.inp
  t0=$(date +%s); echo "[run ] $job  $(date '+%H:%M')"
  ( cd "$W" && ORCA_BIN="$ORCA_BIN" bash run.sh > run.log 2>&1 )
  for f in orca.out step0.out step0a.out step1.out run.log rotate.txt; do [ -f "$W/$f" ] && cp "$W/$f" "$J/"; done
  if [ -f "$J/orca.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$J/orca.out"; then st=OK; rm -rf "$W"; else st=FAILED; tail -3 "$J/run.log"; fi
  echo "[$st] $job  $(( ($(date +%s) - t0) / 60 )) min"
done
