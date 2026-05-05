#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J REGENIE_STEP1
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-2          # 1: Quantitative, 2: Binary

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
activate_env assoc_env

# ---- Parameters ----
TARGET="sv_snv"
step1_r2_cutoff=0.9
step1_maf_cutoff=0.05
step2_r2_cutoff=0.8
step2_maf_cutoff=0.005

PGEN_IN="${OUT_DIR}/${TARGET}_regenie/regenie1"
KEEP_INDV="${OUT_DIR}/sample_selection/twb2_id_train_80.indv"
PHENO="${PHENO_FILE}"
OUTDIR="${OUT_DIR}/${TARGET}_regenie_prs"

QUAN_PHENO=$(get_pheno_list "Quantitative")
BI_PHENO=$(get_pheno_list "Binary")

# MYOPIA_600,MYOPIA_1000 removed due to small sample size.
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_600")
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_1000")
COVAR="AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
CAT_COVAR="SEX"

mkdir -p $OUTDIR $LOG_DIR

# ---- Set TYPE based on array task ID ----
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    TYPE="quan"
    PHENO_LIST=$QUAN_PHENO
    EXTRA_PARAMS="--apply-rint"
else
    TYPE="bi"
    PHENO_LIST=$BI_PHENO
    EXTRA_PARAMS="--bt --minCaseCount 100"
fi

# ---- Log ----
COMBINED_LOG="${LOG_DIR}/02_REGENIE_S1.${TARGET}.${TYPE}.log"
exec > "$COMBINED_LOG" 2>&1

# ---- Run REGENIE Step 1 (all chr merged) ----
echo "Running REGENIE Step 1 | ${TYPE} at $(date)"
regenie \
  --step 1 \
  --keep "${KEEP_INDV}" \
  --pgen "${PGEN_IN}" \
  --phenoFile "${PHENO}" \
  --covarFile "${PHENO}" \
  --phenoColList "${PHENO_LIST}" \
  --covarColList "${COVAR}" \
  --catCovarList "${CAT_COVAR}" \
  $EXTRA_PARAMS \
  --out "${OUTDIR}/regenie_step1_${TYPE}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --bsize 1000 \
  --loocv \
  --lowmem

echo "Step 1 completed | ${TYPE} at $(date)"
echo "Job Array Task $SLURM_ARRAY_TASK_ID ($TYPE) completed at $(date)"
