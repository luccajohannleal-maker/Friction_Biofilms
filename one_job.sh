#!/bin/bash
#SBATCH --job-name=FricBiof_MSc1
#SBATCH --output=/FreeGrow_FrictionBiofilms/log/FricBiof_MSc1%j.txt
#SBATCH --error=/FreeGrow_FrictionBiofilms/log/FricBiof_MSc1%j.txt
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
BASE_DIR="/FreeGrow_FrictionBiofilms/build/Main"

REPEAT_DIR="/data_production/VERTICAL_ORI/lambda${FRIC_COEFF}/repeat${REPEAT_INDEX}/"

# Create the directory for the repeat
mkdir -p "$REPEAT_DIR"

"${BASE_DIR}/main.out" "${REPEAT_DIR}"
