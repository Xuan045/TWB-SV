#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J sum_plt_annot
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
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

# ---- Environment ----
# Environment is initialized in config.sh
activate_env r_env

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/vep
mkdir -p ${logdir}
logfile=${logdir}/06_gwas_summary_plot_annot.log
exec > "$logfile" 2>&1

# ---- Run R script ----
Rscript 06_gwas_summary_plot_annot.r \
    TRUE \
    "${OUT_DIR}/summary_associations/tables"\
    "pleiotropic" 5
