#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J SORT_PHENO
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/../config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/../config.sh"
    elif [ -f "$SLURM_SUBMIT_DIR/scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/scripts/config.sh"
    fi
else
    source "$(dirname "$0")/../config.sh"
fi

# ---- Environment ----
# Environment is initialized in config.sh

# ---- Redirect LOG file ----
mkdir -p "${LOG_DIR}"
logfile="${LOG_DIR}/06_prepare_pheno.log"
exec > "$logfile" 2>&1

echo "Starting phenotype preparation at $(date)"

# ---- Execute Python Script ----
run_in_env data_env python3 06_prepare_pheno.py

echo "Phenotype preparation completed at $(date)"
