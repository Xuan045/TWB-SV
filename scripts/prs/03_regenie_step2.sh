#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J REGENIE_STEP2
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --array=1-44     # 1-22: Quantitative (chr1-22), 23-44: Binary (chr1-22)
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

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
TARGET="sv_snv"
step2_r2_cutoff=0.8
step2_maf_cutoff=0.005

KEEP_INDV="${OUT_DIR}/sample_selection/twb2_id_train_80.indv"
PHENO="${PHENO_FILE}"
JDIR="/work/u4432941/sv"
PFILE_PATH="${JDIR}/${TARGET}_post_imp_v2"
OUTDIR="${OUT_DIR}/${TARGET}_regenie_prs"

QUAN_PHENO=$(get_pheno_list "Quantitative")
BI_PHENO=$(get_pheno_list "Binary")

# MYOPIA_600,MYOPIA_1000 removed due to small sample size.
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_600")
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_1000")
COVAR="AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
CAT_COVAR="SEX"

mkdir -p $OUTDIR $LOG_DIR

# ---- Set TYPE and CHR based on array task ID ----
# Task  1-22: quan, chr1-22
# Task 23-44: bi,   chr1-22
if [ "$SLURM_ARRAY_TASK_ID" -le 22 ]; then
    TYPE="quan"
    CHR=$SLURM_ARRAY_TASK_ID
    PHENO_LIST=$QUAN_PHENO
    EXTRA_PARAMS="--apply-rint"
else
    TYPE="bi"
    CHR=$((SLURM_ARRAY_TASK_ID - 22))
    PHENO_LIST=$BI_PHENO
    EXTRA_PARAMS="--bt --minCaseCount 100"
fi

# ---- Log ----
COMBINED_LOG="${LOG_DIR}/03_REGENIE_S2.${TARGET}.${TYPE}.chr${CHR}.log"
exec > "$COMBINED_LOG" 2>&1

echo "Starting Step 2 | ${TYPE} | chr${CHR} at $(date)"

# ---- pred list from Step 1 output ----
PRED_LIST="${OUTDIR}/regenie_step1_${TYPE}_pred.list"

# ---- Run REGENIE Step 2 ----
regenie \
    --step 2 \
    --pgen ${PFILE_PATH}/chr${CHR}.rsq${step2_r2_cutoff}_maf${step2_maf_cutoff} \
    --keep $KEEP_INDV \
    --chr ${CHR} \
    --bsize 1000 \
    --phenoFile $PHENO \
    --covarFile $PHENO \
    --phenoColList $PHENO_LIST \
    --covarColList $COVAR \
    --catCovarList $CAT_COVAR \
    $EXTRA_PARAMS \
    --pred $PRED_LIST \
    --out ${OUTDIR}/rsq${step2_r2_cutoff}_maf${step2_maf_cutoff}.regenie_step2_${TYPE}_chr${CHR} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --lowmem

echo "Step 2 completed | ${TYPE} | chr${CHR} at $(date)"
echo "Job Array Task $SLURM_ARRAY_TASK_ID (${TYPE} chr${CHR}) completed at $(date)"
