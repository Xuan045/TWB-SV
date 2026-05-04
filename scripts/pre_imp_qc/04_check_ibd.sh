#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J check_ibd
#SBATCH -p ngs372G
#SBATCH -c 56
#SBATCH --mem=350g
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# ---- Source Configuration ----
source "$(dirname "$0")/../config.sh"

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

KING="/staging/biology/u4432941/apps/king"

# ---- Redirect LOG file ----
logfile=${LOG_DIR}/04_check_ibd.log
exec > "$logfile" 2>&1

TWB1=${OUT_DIR}/pre_imp_qc_v1/twb1_preqc_sexcheck_hetcheck
TWB2=${OUT_DIR}/pre_imp_qc_v2/twb2_preqc_sexcheck_hetcheck
ibddir=${OUT_DIR}/check_ibd
mkdir -p $ibddir

# ── Within-batch kinship ──────────────────────────────────────────
for batch in 1 2; do
    if [ $batch -eq 1 ]; then
        BFILE=$TWB1
    else
        BFILE=$TWB2
    fi
    PARA=${ibddir}/twb${batch}

    # LD pruning
    plink --bfile $BFILE \
        --autosome \
        --maf 0.05 \
        --snps-only just-acgt \
        --indep-pairwise 200 100 0.1 \
        --out ${PARA}.ldpr

    # Within-batch kinship estimation
    plink2 --bfile $BFILE \
        --extract ${PARA}.ldpr.prune.in \
        --make-king-table \
        --king-table-filter 0.0884 \
        --out ${PARA}.ldpr.2degree


    plink --bfile $BFILE \
        --extract ${PARA}.ldpr.prune.in \
        --make-bed \
        --out ${PARA}.ldpr
    $KING -b ${PARA}.ldpr.bed \
        --related \
        --degree 2 \
        --prefix ${PARA}.king

done

# ── Between-batch kinship ──────────────────────────────────────────
PARA_MERGE=${ibddir}/twb_merged

# Find common SNPs
awk 'NR==FNR{a[$2];next} ($2 in a){print $2}' \
    ${TWB1}.bim ${TWB2}.bim > ${PARA_MERGE}_comm.snplist

# Extract common SNPs
plink --bfile $TWB1 --extract ${PARA_MERGE}_comm.snplist --make-bed --out ${TWB1}_tmp
plink --bfile $TWB2 --extract ${PARA_MERGE}_comm.snplist --make-bed --out ${TWB2}_tmp

# Merge
plink --bfile ${TWB1}_tmp \
    --bmerge ${TWB2}_tmp \
    --make-bed \
    --out ${PARA_MERGE}

# LD pruning on merged
plink --bfile ${PARA_MERGE} \
    --autosome \
    --maf 0.05 \
    --snps-only just-acgt \
    --indep-pairwise 200 100 0.1 \
    --out ${PARA_MERGE}.ldpr

# Across-batch kinship estimation
plink2 --bfile ${PARA_MERGE} \
    --extract ${PARA_MERGE}.ldpr.prune.in \
    --make-king-table \
    --king-table-filter 0.0884 \
    --out ${PARA_MERGE}.ldpr.2degree

# KING
plink --bfile ${PARA_MERGE} \
    --extract ${PARA_MERGE}.ldpr.prune.in \
    --make-bed \
    --out ${PARA_MERGE}.ldpr
$KING -b ${PARA_MERGE}.ldpr.bed \
    --related \
    --degree 2 \
    --prefix ${PARA_MERGE}.king

rm ${TWB1}_tmp* ${TWB2}_tmp*
# Remove temporary ldpr bfiles
rm ${PARA_MERGE}.ldpr.{bim,bed,fam} ${ibddir}/twb1.ldpr.{bim,bed,fam} ${ibddir}/twb2.ldpr.{bim,bed,fam}