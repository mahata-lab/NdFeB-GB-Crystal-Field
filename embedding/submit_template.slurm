#!/bin/bash
# SLURM template: copy into a job folder (or loop over folders) and set ORCA_DIR and the module line for your cluster.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --job-name=orca_4f
#SBATCH --output=run_%j.log
module load openmpi/4.1.8                     # OpenMPI version matching the ORCA 6.1.1 build
ORCA_DIR=/path/to/orca_6_1_1                  # <-- edit
export ORCA_BIN="$ORCA_DIR/orca" PATH="$ORCA_DIR:$PATH" LD_LIBRARY_PATH="$ORCA_DIR:${LD_LIBRARY_PATH:-}" OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
bash run.sh
