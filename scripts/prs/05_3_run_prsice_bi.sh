#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J PRSICE_BI
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --array=1-46          # To run all binary traits, change to --array=1-46
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

# ---- Parameters (passed via --export) ----
TARGET=${TARGET:?Error: TARGET not set. Pass via --export=ALL,TARGET=sv_snv}

BI_PHENO=$(get_pheno_list "Binary")

# Exclude MYOPIA_600 and MYOPIA_1000 since they were removed in step 1/2 due to small sample size
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_600")
BI_PHENO=$(exclude_pheno "$BI_PHENO" "MYOPIA_1000")

IFS=',' read -ra PHENO_LIST <<< "$BI_PHENO"
PHENO_NAME=${PHENO_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

PRSICE="${PRSICE_BIN}"
regenie_dir="${OUT_DIR}/${TARGET}.rsq0.8_maf0.005.regenie_results_prs"
outdir="${OUT_DIR}/prsice/${TARGET}.bi"
LOG_DIR="${PRSICE_LOG_DIR}"

MERGED="${OUT_DIR}/prsice/rsq0.8_maf0.005.merged"
# PHENO_FILE is defined in config.sh
TEST_INDV="${OUT_DIR}/sample_selection/twb2_id_test_20.indv"
SUM_STATS="${regenie_dir}/${PHENO_NAME}.regenie"

mkdir -p $LOG_DIR $outdir

# ---- Log ----
COMBINED_LOG="${LOG_DIR}/05_prsice.${TARGET}.bi.${PHENO_NAME}.log"
exec > "$COMBINED_LOG" 2>&1

echo "Starting PRSice-2 | ${TARGET} | bi | ${PHENO_NAME} at $(date)"

# ---- Reformat sumstats ----
activate_env data_env
python3 05_1_reformat_prsice_table.py \
    --sumstats $SUM_STATS \
    --out ${outdir}/${PHENO_NAME}_sumstats_for_prsice.txt

# ---- Run PRSice-2 ----
activate_env assoc_env

$PRSICE \
    --base ${outdir}/${PHENO_NAME}_sumstats_for_prsice.txt \
    --snp SNP --a1 A1 --a2 A2 --beta --stat BETA --se SE --pvalue P \
    --target ${MERGED} \
    --keep $TEST_INDV \
    --pheno $PHENO_FILE \
    --pheno-col $PHENO_NAME \
    --cov $PHENO_FILE \
    --cov-col AGE,SEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
    --cov-factor SEX \
    --binary-target T \
    --prevalence 0.1 \
    --thread $SLURM_CPUS_PER_TASK \
    --out ${outdir}/${PHENO_NAME}

echo "PRSice-2 bi completed | ${PHENO_NAME} at $(date)"
