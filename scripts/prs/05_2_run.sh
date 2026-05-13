#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J PRS_CV
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH --array=10-12
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

# ---- Parameters ----
QUAN_PHENO=$(get_pheno_list "Quantitative")
IFS=',' read -ra PHENO_LIST <<< "$QUAN_PHENO"
PHENO_NAME=${PHENO_LIST[$((SLURM_ARRAY_TASK_ID - 1))]}

LOG_DIR="${PRSICE_LOG_DIR}"
mkdir -p "$LOG_DIR"

# ---- Log: all stdout/stderr from this point go to a single log ----
COMBINED_LOG="${LOG_DIR}/05_prs_repeats.${PHENO_NAME}.log"
exec > "$COMBINED_LOG" 2>&1

echo "========================================================"
echo " PRS Repeats Pipeline | ${PHENO_NAME} | $(date)"
echo "========================================================"

# ---- Step 1: PRSice-2 for sv_snv ----
echo ""
echo "[Step 1/3] PRSice-2 | sv_snv | ${PHENO_NAME}"
bash "${SLURM_SUBMIT_DIR}/05_2_prsice_quan.sh" "$PHENO_NAME" "sv_snv"
if [ $? -ne 0 ]; then
    echo "ERROR: PRSice-2 failed for sv_snv | ${PHENO_NAME}" >&2
    exit 1
fi

# ---- Step 2: PRSice-2 for snv ----
echo ""
echo "[Step 2/3] PRSice-2 | snv | ${PHENO_NAME}"
bash "${SLURM_SUBMIT_DIR}/05_2_prsice_quan.sh" "$PHENO_NAME" "snv"
if [ $? -ne 0 ]; then
    echo "ERROR: PRSice-2 failed for snv | ${PHENO_NAME}" >&2
    exit 1
fi

# ---- Step 3: 100-repeat CV and comparison ----
echo ""
echo "[Step 3/3] 100-repeat CV | ${PHENO_NAME}"
activate_env r_env
Rscript "${SLURM_SUBMIT_DIR}/05_2_prsice_repeats.r" \
    "$PHENO_NAME" \
    "sv_snv" \
    "snv" \
    100 \
    5000

if [ $? -ne 0 ]; then
    echo "ERROR: CV script failed | ${PHENO_NAME}" >&2
    exit 1
fi

echo ""
echo "========================================================"
echo " Done | ${PHENO_NAME} | $(date)"
echo "========================================================"
