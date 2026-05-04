#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J FILTER_SIG       # Job name
#SBATCH -p ngs53G           # Partition Name 
#SBATCH -c 8               
#SBATCH --mem=53G           # memory used
#SBATCH -o out.log
#SBATCH -e /dev/null

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

target_panel="sv_snv" # sv_snv or snv
FILEPATH="${OUT_DIR}/${target_panel}.rsq0.8_maf0.005.regenie_results"
OUTPATH="${FILEPATH}/filtered_strict"
# LOG10P_THRESHOLD=7.30103
# 5*10^(-8)/86 traits
LOG10P_THRESHOLD=9.235528

mkdir -p "$OUTPATH"


awk -F, 'NR>1 {print $1}' "$PHENO_META" | while read -r pheno; do
  echo "Processing phenotype: $pheno"

  INPUT="${FILEPATH}/${pheno}.regenie"
  OUTPUT="${OUTPATH}/${pheno}.filtered.regenie"

  awk -v threshold="$LOG10P_THRESHOLD" 'BEGIN { FS=OFS="\t" }
  NR==1 { print; next }
  {
    if ($12 > threshold)
      print
  }' "$INPUT" > "$OUTPUT"

done
