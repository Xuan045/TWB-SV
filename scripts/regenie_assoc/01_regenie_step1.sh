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
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# ---- Parameters ----
TARGET="snv"
step1_r2_cutoff=0.9
step1_maf_cutoff=0.05
step2_r2_cutoff=0.8
step2_maf_cutoff=0.005

WKDIR=${PROJECT_DIR}
JDIR=${WORK_DIR:-/work/u4432941/sv}/${TARGET}_post_imp_v2
OUT=${OUT_DIR}/${TARGET}_regenie
LOG_DIR="${LOG_DIR:-${WKDIR}/logs}/regenie"
PHENO="$PHENO_FILE"

QUAN_PHENO=$(get_pheno_list "Quantitative")
BI_PHENO=$(get_pheno_list "Binary")

COVAR="AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
CAT_COVAR="SEX"

mkdir -p $OUT $LOG_DIR

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
COMBINED_LOG="${LOG_DIR}/REGENIE_S1.${TARGET}.${TYPE}.log"
exec > "$COMBINED_LOG" 2>&1

# ---- Prepare shared files (only task 1) ----
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "Task 1: Preparing shared variant list and merge list..."

    rm -f ${OUT}/regenie1.varlist.tmp
    for chr in {1..22}; do
        cat "${JDIR}/${TARGET}.rsq${step1_r2_cutoff}_maf${step1_maf_cutoff}_snv/chr${chr}.rsq${step1_r2_cutoff}_maf${step1_maf_cutoff}_snv.ldpr.geno0.05.prune.in" >> ${OUT}/regenie1.varlist.tmp
    done
    mv ${OUT}/regenie1.varlist.tmp ${OUT}/regenie1.varlist

    rm -f ${OUT}/pgen.merge_list.tmp
    for chr in {1..22}; do
        echo "${JDIR}/chr${chr}.rsq${step2_r2_cutoff}_maf${step2_maf_cutoff}" >> ${OUT}/pgen.merge_list.tmp
    done
    mv ${OUT}/pgen.merge_list.tmp ${OUT}/pgen.merge_list.txt

    echo "Task 1: Shared files are ready."
else
    echo "Task 2: Waiting for shared files..."
    while [ ! -f "${OUT}/pgen.merge_list.txt" ] || [ ! -f "${OUT}/regenie1.varlist" ]; do
        sleep 5
    done
    echo "Task 2: Shared files detected. Proceeding..."
fi

sleep 10s

# ---- Generate merged PGEN across all chr (only task 1 runs plink2) ----
if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    if [ ! -f "${OUT}/regenie1.pgen" ]; then
        echo "Task 1: Generating merged PGEN for all chromosomes..."
        plink2 \
          --pmerge-list "${OUT}/pgen.merge_list.txt" \
          --extract "${OUT}/regenie1.varlist" \
          --maf 0.05 \
          --geno 0.05 \
          --make-pgen \
          --out "${OUT}/regenie1"
    fi
else
    echo "Task 2: Waiting for merged PGEN..."
    while true; do
        if [ -f "${OUT}/regenie1.log" ] && grep -q "End time" "${OUT}/regenie1.log"; then
            break
        fi
        sleep 10
    done
    sleep 5
fi

# ---- Run REGENIE Step 1 (all chr merged) ----
echo "Running REGENIE Step 1 | ${TYPE} at $(date)"
regenie \
  --step 1 \
  --pgen "${OUT}/regenie1" \
  --phenoFile "${PHENO}" \
  --covarFile "${PHENO}" \
  --phenoColList "${PHENO_LIST}" \
  --covarColList "${COVAR}" \
  --catCovarList "${CAT_COVAR}" \
  $EXTRA_PARAMS \
  --out "${OUT}/regenie_step1_${TYPE}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --bsize 1000 \
  --loocv \
  --lowmem

echo "Step 1 completed | ${TYPE} at $(date)"
echo "Job Array Task $SLURM_ARRAY_TASK_ID ($TYPE) completed at $(date)"
