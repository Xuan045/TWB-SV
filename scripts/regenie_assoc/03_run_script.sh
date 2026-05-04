#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J REGENIE_merge  # Job name
#SBATCH -p ngs53G           # Partition Name 
#SBATCH -c 8               
#SBATCH --mem=53G           # memory used
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-86

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh

PHENO=$(awk -F, 'NR>1 {print $1}' "$PHENO_META" | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Processing phenotype: $PHENO"

# -----------------------------
# Run combine_results.py
# -----------------------------
activate_env data_env

python 03_1_combine_results.py --phenotype "$PHENO" --panel "sv_snv"
python 03_1_combine_results.py --phenotype "$PHENO" --panel "snv"

# -----------------------------
# Run 03_2_regenie_plot.r
# -----------------------------
activate_env r_env

# The second argument is the boolean to indicate whether to plot the results restricted to SV
# SNV+SV
Rscript 03_2_regenie_plot.r $PHENO sv_snv FALSE
Rscript 03_2_regenie_plot.r $PHENO sv_snv TRUE

# SNV only
Rscript 03_2_regenie_plot.r $PHENO snv FALSE
