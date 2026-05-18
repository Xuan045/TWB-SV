#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J VEP_SNV
#SBATCH -p ngs186G
#SBATCH -c 28
#SBATCH --mem=175g
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-22

# ---- Source Configuration ----
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    if [ -f "$SLURM_SUBMIT_DIR/../../scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/../../scripts/config.sh"
    elif [ -f "$SLURM_SUBMIT_DIR/scripts/config.sh" ]; then
        source "$SLURM_SUBMIT_DIR/scripts/config.sh"
    fi
else
    source "$(dirname "$0")/../../scripts/config.sh"
fi

# ---- Redirect LOG file ----
logdir=${LOG_DIR}/vep
mkdir -p ${logdir}
logfile=${logdir}/03_vep_snv_${SLURM_ARRAY_TASK_ID}.log
exec > "$logfile" 2>&1

# ---- Environment ----
# Environment is initialized in config.sh
activate_env assoc_env

# target imputation panel should be either "sv_snv" or "snv"
TARGET="sv_snv"
CHR=${SLURM_ARRAY_TASK_ID}
# CHR=1

WKDIR=${OUT_DIR}
PFILE_PATH=${WORK_DIR}/${TARGET}_post_imp_v2
OUT=$WKDIR/vep/snv
SIG_LIST=$WKDIR/summary_associations/varlists/quan/chr${CHR}_snv.varlist
SAMPLE_ID=rsq0.8_maf0.005

mkdir -p $OUT

# -------------------------------------------------
# (Filter to significant variants and) make VCF file
# -------------------------------------------------
plink2 \
    --pfile ${PFILE_PATH}/chr${CHR}.${SAMPLE_ID} \
    --extract $SIG_LIST \
    --make-pgen pvar-cols=vcfheader,qual,filter,info \
    --output-chr chr26 \
    --out ${OUT}/temp.chr${CHR}

# Rename the pvar file to vcf
mv ${OUT}/temp.chr${CHR}.pvar ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf
bcftools view -Oz -o ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf.gz ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf
tabix -p vcf ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf.gz

rm ${OUT}/temp.chr${CHR}.pgen ${OUT}/temp.chr${CHR}.psam ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf

INPUT_VCF_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.vcf.gz
OUTPUT_VCF_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.vep.vcf.gz
OUTPUT_TSV_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.tsv

activate_env $VEP_ENV_V115

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started" >> ${logfile}

custom_twb="${TWB_DRAGEN},TWB1490_WGS,vcf,exact,0,AF"

vep --cache --offline \
    -i $INPUT_VCF_PATH \
    --format vcf \
    --fork 4 \
    --check_existing \
    --force_overwrite \
    --dir_cache $VEP_CACHE_DIR \
    --assembly GRCh38 \
    --overlaps \
    --pick \
    --plugin GWAS,file=$VEP_GWAS \
    --custom ${custom_twb} \
    --no_stats \
    --fasta $FASTA_REF \
    --vcf \
    -o $OUTPUT_VCF_PATH

# generate tsv
echo -e "CHROM\tPOS\tID\tREF\tALT\t$(bcftools +split-vep -l $OUTPUT_VCF_PATH | cut -f2 | tr '\n' '\t' | sed 's/\t$//')" > $OUTPUT_TSV_PATH
bcftools +split-vep -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\n' -d -A tab $OUTPUT_VCF_PATH >> $OUTPUT_TSV_PATH

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished" >> ${logfile}
