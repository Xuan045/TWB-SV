#!/usr/bin/bash
# ==============================================================================
# Run PRSice-2 for a single phenotype and target panel.
#
# Usage:
#     bash 05_1_prsice_run.sh <phenotype> <target>
#
# Outputs:
#     <OUT_DIR>/prsice/<target>.quan/<phenotype>.all_score
#     <OUT_DIR>/prsice/<target>.quan/<phenotype>.log
# ==============================================================================

PHENO_NAME=${1:?Error: phenotype not provided}
TARGET=${2:?Error: target not provided}

# ---- Source Configuration ----
# Must be sourced explicitly since this script is called as a subprocess
SCRIPT_DIR="$(dirname "$0")"
if [ -f "${SCRIPT_DIR}/../config.sh" ]; then
    source "${SCRIPT_DIR}/../config.sh"
elif [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "ERROR: config.sh not found relative to ${SCRIPT_DIR}"
    exit 1
fi

PRSICE="${PRSICE_BIN}"
regenie_dir="${OUT_DIR}/${TARGET}.rsq0.8_maf0.005.regenie_results_prs"
outdir="${OUT_DIR}/prsice/${TARGET}.quan"
MERGED="${OUT_DIR}/prsice/test.${TARGET}.rsq0.8_maf0.005.merged"
TEST_INDV="${OUT_DIR}/sample_selection/twb2_id_test_20.indv"
SUM_STATS="${regenie_dir}/${PHENO_NAME}.regenie"
REFORMATTED="${outdir}/${PHENO_NAME}_sumstats_for_prsice.txt"

mkdir -p "$outdir"

# ---- Reformat sumstats ----
echo "  Reformatting sumstats | ${TARGET} | ${PHENO_NAME}"
activate_env data_env
python3 "$(dirname "$0")/05_1_reformat_prsice_table.py" \
    --sumstats "$SUM_STATS" \
    --out "$REFORMATTED"

if [ $? -ne 0 ]; then
    echo "ERROR: Reformat failed for ${TARGET} | ${PHENO_NAME}"
    exit 1
fi

# ---- Run PRSice-2 ----
# --all-score : output PRS at every p-value threshold
# --no-regress: skip PRSice-2 internal regression; evaluation handled in R
echo "  Running PRSice-2 | ${TARGET} | ${PHENO_NAME}"
activate_env r_env

$PRSICE \
    --base "$REFORMATTED" \
    --chr CHR --bp POS --snp SNP \
    --a1 A1 --a2 A2 --beta --stat BETA --pvalue P \
    --target "$MERGED" \
    --keep "$TEST_INDV" \
    --pheno "$PHENO_FILE" \
    --pheno-col "$PHENO_NAME" \
    --cov "$PHENO_FILE" \
    --cov-col AGE,SEX,@PC[1-10] \
    --cov-factor SEX \
    --binary-target F \
    --all-score \
    --no-regress \
    --thread 14 \
    --out "${outdir}/${PHENO_NAME}"

if [ $? -ne 0 ]; then
    echo "ERROR: PRSice-2 failed for ${TARGET} | ${PHENO_NAME}"
    exit 1
fi

echo "  PRSice-2 done | ${TARGET} | ${PHENO_NAME}"
