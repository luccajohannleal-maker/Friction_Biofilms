#!/bin/bash
#SBATCH --job-name=FricBiof_MSc
#SBATCH --output=./log/FricBiof_MSc%j.txt
#SBATCH --error=./log/FricBiof_MSc%j.txt
#SBATCH --partition=long
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB

export OMP_NUM_THREADS="16"
ulimit -c unlimited # this is for debugging if needed

# Get the repeat index from the argument
REPEAT_INDEX=$1
FRIC_COEFF=8


# Define directories
BASE_DIR="./build/Main"

REPEAT_DIR="data_production/Lambda${FRIC_COEFF}/repeat${REPEAT_INDEX}/"

# Create the directory for the repeat
mkdir -p "$REPEAT_DIR"

"${BASE_DIR}/main.out" "${REPEAT_DIR}"
