#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J PGEN_MERGE
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
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

pgen_dir="/work/u4432941/sv/${TARGET}_post_imp_v2"
outdir="${OUT_DIR}/prsice"
LOG_DIR="${PRSICE_LOG_DIR}"
TEST_INDV="${OUT_DIR}/sample_selection/twb2_id_test_20.indv"

mkdir -p $LOG_DIR $outdir

# ---- Log ----
COMBINED_LOG="${LOG_DIR}/05_pgen_merge.${TARGET}.log"
exec > "$COMBINED_LOG" 2>&1

echo "Starting merge | ${TARGET} at $(date)"

# ---- Merge pgen and combine to bim ----
MERGED=${outdir}/test.${TARGET}.rsq0.8_maf0.005.merged

for chr in {1..22}; do
    echo "${pgen_dir}/chr${chr}.rsq0.8_maf0.005"
done > ${outdir}/${TARGET}_merge_list.txt

plink2 \
    --pmerge-list ${outdir}/${TARGET}_merge_list.txt \
    --keep $TEST_INDV \
    --make-bed \
    --out ${MERGED}

echo "Merge completed at $(date)"
