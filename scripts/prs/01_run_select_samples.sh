#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J RUN_SELECT
#SBATCH -p ngs186G
#SBATCH -c 28
#SBATCH --mem=175g
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

# ---- Log ----
mkdir -p "$LOG_DIR"

LOG="${LOG_DIR}/01_run_select.log"
exec > "$LOG" 2>&1

# ---- Tool dependencies ----
ml biology bcftools/1.13

# Define local variables based on config
OUTDIR="${OUT_DIR}/sample_selection"
ID_OUT="${OUTDIR}/valid_twb2_id.txt"
OUT_PREFIX="${OUTDIR}/twb2_id"
mkdir -p "$OUTDIR"

# ---- Extract and split sample IDs ----
run_in_env data_env \
  python3 01_1_select_samples.py \
    "$LAB_INFO" \
    "$ID_OUT" \
    "$BCF_PHASED_V2" \
    "$OUT_PREFIX"

# ---- Evaluate kinship coefficients ----
activate_env assoc_env

BFILE="$BFILE_PCA"
TRAIN_INDV="${OUT_PREFIX}_train_80.indv"
TEST_INDV="${OUT_PREFIX}_test_20.indv"
OUT_KIN="${OUTDIR}/kinship"
mkdir -p "$OUT_KIN"

# Link v2 relatedness files to work_dir
ln -sf "$KIN0_FILE" "$OUT_KIN/twb2.ldpr.2degree.kin0"

plink2 \
  --bfile $BFILE \
  --keep $TEST_INDV \
  --autosome \
  --maf 0.05 \
  --snps-only just-acgt \
  --indep-pairwise 200 100 0.1 \
  --out ${OUT_KIN}/ldpr

# Evaluate within test set relatedness
plink2 \
  --bfile $BFILE \
  --keep $TEST_INDV \
  --extract ${OUT_KIN}/ldpr.prune.in \
  --king-cutoff 0.0884 \
  --out $OUT_KIN/twb2_test.2degree

# ---- Final pruning related samples ----
plink2 \
  --bfile $BFILE \
  --remove $OUT_KIN/twb2_test.2degree.king.cutoff.out.id \
  --autosome \
  --maf 0.05 \
  --snps-only just-acgt \
  --indep-pairwise 200 100 0.1 \
  --out ${OUT_KIN}/twb2_rel_test_removed.ldpr

# Evaluate within test set relatedness
plink2 \
  --bfile $BFILE \
  --remove $OUT_KIN/twb2_test.2degree.king.cutoff.out.id \
  --extract ${OUT_KIN}/twb2_rel_test_removed.ldpr.prune.in \
  --make-king-table \
  --king-table-filter 0.0884 \
  --out $OUT_KIN/twb2_rel_test_removed.2degree

# ---- Generate final remove list ----
run_in_env data_env \
  python3 01_2_final_remove_list.py \
    --pairs $OUT_KIN/twb2_rel_test_removed.2degree.txt \
    --test ${OUT_PREFIX}_test_20.txt \
    --train ${OUT_PREFIX}_train_80.txt \
    --output $OUTDIR/twb2_final
