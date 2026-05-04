#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J REGENIE_summary_plot  # Job name
#SBATCH -p ngs53G           # Partition Name 
#SBATCH -c 8               
#SBATCH --mem=53G           # memory used


# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env r_env

# ---- Set OUTPUT_DIR for R script ----
# Uses OUT_DIR from config.sh
export OUTPUT_DIR="${OUT_DIR}"

# -----------------------------
# Run 05_summary_plot.r
# -----------------------------
Rscript 05_summary_plot.r

# Filter to SV
Rscript 05_summary_plot.r TRUE