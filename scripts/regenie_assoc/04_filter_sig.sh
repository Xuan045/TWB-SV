#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J FILTER_SIG       # Job name
#SBATCH -p ngs53G           # Partition Name 
#SBATCH -c 8               
#SBATCH --mem=53G           # memory used
#SBATCH -o out.log
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

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/regenie
mkdir -p ${logdir}
logfile=${logdir}/04_filter_sig.log
exec > "$logfile" 2>&1

target_panel="sv_snv" # sv_snv or snv
FILEPATH="${OUT_DIR}/${target_panel}.rsq0.8_maf0.005.regenie_results"

# Define the two output directories
OUTPATH_STRICT="${FILEPATH}/filtered_strict"
OUTPATH_5E8="${FILEPATH}/filtered_5e-8"

mkdir -p "$OUTPATH_STRICT" "$OUTPATH_5E8"

# Thresholds
# 5*10^(-8)
LOG10P_THRESHOLD=7.30103
# 5*10^(-8)/86 traits
LOG10P_THRESHOLD_STRICT=9.235528

awk -F, 'NR>1 {print $1}' "$PHENO_META" | while read -r pheno; do
  echo "Processing phenotype: $pheno"

  INPUT="${FILEPATH}/${pheno}.regenie"
  
  # Check if input exists to avoid errors
  if [[ ! -f "$INPUT" ]]; then
    echo "Warning: $INPUT not found. Skipping..."
    continue
  fi

  # Use awk to split the stream into two different files based on logic
  awk -v t_5e8="$LOG10P_THRESHOLD" \
      -v t_strict="$LOG10P_THRESHOLD_STRICT" \
      -v out_5e8="${OUTPATH_5E8}/${pheno}.filtered.regenie" \
      -v out_strict="${OUTPATH_STRICT}/${pheno}.filtered.regenie" \
      'BEGIN { FS=OFS="\t" }
      NR==1 { 
          print > out_5e8; 
          print > out_strict; 
          next 
      }
      {
          # If it passes the lower threshold (5e-8), print to the 5e-8 folder
          if ($12 > t_5e8) {
              print > out_5e8
          }
          # If it also passes the higher threshold (strict), print to the strict folder
          if ($12 > t_strict) {
              print > out_strict
          }
      }' "$INPUT"

done
