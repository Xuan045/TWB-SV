#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J merge
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
activate_env data_env

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/vep
mkdir -p ${logdir}
logfile=${logdir}/04_merge.log
exec > "$logfile" 2>&1

# ---- Run Python script ----
python3 04_merge_vep_annotations.py
python3 05_filter_gwas_variants.py
