#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J PRSICE
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --array=1-38          # To run all qunatitative traits, change to --array=1-38
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh

# ---- Parameters (passed via --export) ----
TARGET=${TARGET:?Error: TARGET not set. Pass via --export=ALL,TARGET=sv_snv}

QUAN_PHENO=$(get_pheno_list "Quantitative")

IFS=',' read -ra PHENO_LIST <<< "$QUAN_PHENO"
PHENO_NAME=${PHENO_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

PRSICE="${PRSICE_BIN}"
pgen_dir="/work/u4432941/sv/${TARGET}_post_imp_v2"
regenie_dir="${OUT_DIR}/${TARGET}.rsq0.8_maf0.005.regenie_results_prs"
outdir="${OUT_DIR}/prsice/${TARGET}.quan"
LOG_DIR="${PRSICE_LOG_DIR}"

MERGED="${OUT_DIR}/prsice/test.${TARGET}.rsq0.8_maf0.005.merged"
# PHENO_FILE is defined in config.sh
TEST_INDV="${OUT_DIR}/sample_selection/twb2_id_test_20.indv"
TRAIN_INDV="${OUT_DIR}/sample_selection/twb2_id_train_80.indv"
SUM_STATS="${regenie_dir}/${PHENO_NAME}.regenie"

mkdir -p $LOG_DIR $outdir

# ---- Log ----
COMBINED_LOG="${LOG_DIR}/05_prsice.${TARGET}.quan.${PHENO_NAME}.log"
exec > "$COMBINED_LOG" 2>&1

echo "Starting PRSice-2 | ${TARGET} | quan | ${PHENO_NAME} at $(date)"

# ---- Reformat sumstats for PRSice-2 ----
# PRSice-2 wants: SNP, A1, A2, BETA, SE, P
activate_env data_env
python3 05_1_reformat_prsice_table.py \
    --sumstats $SUM_STATS \
    --out ${outdir}/${PHENO_NAME}_sumstats_for_prsice.txt

# ---- Run PRSice-2 ----
activate_env r_env

$PRSICE \
    --base ${outdir}/${PHENO_NAME}_sumstats_for_prsice.txt \
    --chr CHR --bp POS --snp SNP  \
    --a1 A1 --a2 A2 --beta --stat BETA --pvalue P \
    --target ${MERGED} \
    --keep $TEST_INDV \
    --pheno $PHENO_FILE \
    --pheno-col $PHENO_NAME \
    --cov $PHENO_FILE \
    --cov-col AGE,SEX,@PC[1-10] \
    --cov-factor SEX \
    --binary-target F \
    --perm 10000 \
    --thread $SLURM_CPUS_PER_TASK \
    --out ${outdir}/${PHENO_NAME}

echo "PRSice-2 quan completed | ${PHENO_NAME} at $(date)"
