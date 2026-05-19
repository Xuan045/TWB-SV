#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J REGENIE_summary_plot
#SBATCH -p ngs53G
#SBATCH -c 8
#SBATCH --mem=53G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/../../scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/../../scripts/config.sh"
    elif [ -f "$SLURM_SUBMIT_DIR/scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/scripts/config.sh"
    fi
else
    source "$(dirname "$0")/../../scripts/config.sh"
fi

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/post_assoc
mkdir -p ${logdir}
logfile=${logdir}/summary_plot.log
exec > "$logfile" 2>&1

# ---- Environment ----
# Environment is initialized in config.sh
activate_env r_env

# ---- Set OUTPUT_DIR for R script ----
export OUTPUT_DIR="${OUT_DIR}"

# -----------------------------
# Run summary_plot.r
# -----------------------------
Rscript summary_plot.r

# Filter to SV
Rscript summary_plot.r TRUE