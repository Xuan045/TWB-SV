#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J RUN_COMPARE
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

LOG_DIR="${PRSICE_LOG_DIR}"
mkdir -p "$LOG_DIR"

# ---- Log: all stdout/stderr from this point go to a single log ----
COMBINED_LOG="${LOG_DIR}/06_plotting.log"
exec > "$COMBINED_LOG" 2>&1

activate_env r_env
Rscript 06_plot_comparison.r
