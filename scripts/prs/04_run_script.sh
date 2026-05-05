#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J REGENIE_merge  # Job name
#SBATCH -p ngs53G           # Partition Name 
#SBATCH -c 8               
#SBATCH --mem=53G           # memory used
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-84

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

# PHENO_LIST is defined in config.sh
PHENO=$(sed -n "$((SLURM_ARRAY_TASK_ID))p" $PHENO_LIST)

echo "Processing phenotype: $PHENO"

# -----------------------------
# Run combine_results.py
# -----------------------------
activate_env data_env

python 04_1_combine_results.py --phenotype "$PHENO" --panel "sv_snv"
python 04_1_combine_results.py --phenotype "$PHENO" --panel "snv"

# -----------------------------
# Run 04_2_regenie_plot.r
# -----------------------------
activate_env r_env

# The second argument is the boolean to indicate whether to plot the results restricted to SV
# SNV+SV
Rscript 04_2_regenie_plot.r $PHENO sv_snv FALSE
Rscript 04_2_regenie_plot.r $PHENO sv_snv TRUE

# SNV only
Rscript 04_2_regenie_plot.r $PHENO snv FALSE
