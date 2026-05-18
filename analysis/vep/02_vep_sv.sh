#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J VEP_SV
#SBATCH -p ngs92G
#SBATCH -c 14
#SBATCH --mem=92g
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-22             # chr1-22

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
logfile=${logdir}/02_vep_sv_${SLURM_ARRAY_TASK_ID}.log
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
OUT=$WKDIR/vep/sv
SIG_LIST=$WKDIR/summary_associations/varlists/quan/chr${CHR}_sv.varlist
SAMPLE_ID=rsq0.8_maf0.005

mkdir -p $OUT

# -------------------------------------------------
# (Filter to significant variants and) make VCF file
# -------------------------------------------------
plink2 \
    --pfile ${PFILE_PATH}/chr${CHR}.rsq0.8_maf0.005 \
    --extract $SIG_LIST \
    --make-pgen pvar-cols=vcfheader,qual,filter,info \
    --output-chr chr26 \
    --out ${OUT}/temp.chr${CHR}

# -----------------------------------------------
# Reformat the VCF
# The script will fill in SVLEN and SVTYPE
# based on the sv_info.tsv
# -----------------------------------------------
activate_env data_env
python3 02_add_svlen_svtype.py sv_info.tsv ${OUT}/temp.chr${CHR}.pvar ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf ${logfile}
bgzip -f ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf
tabix -p vcf ${OUT}/chr${CHR}.${SAMPLE_ID}.vcf.gz

rm ${OUT}/temp.chr${CHR}.p*

INPUT_VCF_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.vcf.gz
OUTPUT_VCF_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.vep.vcf.gz
OUTPUT_TSV_PATH=${OUT}/chr${CHR}.${SAMPLE_ID}.tsv

activate_env $VEP_ENV_V112

echo "$(date '+%Y-%m-%d %H:%M:%S') Job started" >> ${logfile}

custom_del="file=${TWB_DEL},short_name=TWB1480_DEL,format=vcf,fields=AF,type=overlap,distance=200"
custom_ins="file=${TWB_INS},short_name=TWB1480_INS,format=vcf,fields=AF,type=overlap,distance=200"

field="FILTER%AF%GRPMAX_AF%AF_eas%AF_non_neuro%AF_non_neuro_eas%AF_controls_and_biobanks%AF_controls_and_biobanks_eas"
gnomad_del="file=${GNOMAD_DEL},short_name=gnomAD_v4.1_overlap100,format=vcf,fields=${field},type=overlap,overlap_cutoff=100,reciprocal=1"
gnomad_ins_1="file=${GNOMAD_INS},short_name=gnomAD_v4.1_dist200,format=vcf,fields=${field},type=overlap,distance=200"
gnomad_ins_2="file=${GNOMAD_INS},short_name=gnomAD_v4.1_overlap50,format=vcf,fields=${field},type=overlap,overlap_cutoff=50,reciprocal=1"
gnomad_ins_3="file=${GNOMAD_INS},short_name=gnomAD_v4.1_exact,format=vcf,fields=${field},type=exact"

clinvar_del="file=${CLINVAR_DEL},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNREVSTAT%MC%CLNDN,type=exact"
clinvar_ins="file=${CLINVAR_INS},short_name=ClinVar,format=vcf,fields=CLNSIG%CLNREVSTAT%MC%CLNDN,type=exact"
clinvar_del_p="file=${CLINVAR_DEL_P},short_name=ClinVar_p_lp,format=vcf,fields=CLNSIG%CLNREVSTAT%MC%CLNDN,type=within"
clinvar_del_b="file=${CLINVAR_DEL_B},short_name=ClinVar_b_lb,format=vcf,fields=CLNSIG%CLNREVSTAT%MC%CLNDN,type=surrounding"

func_impact="file=${ENCODE},short_name=ENCODE_cCREs,format=bed,type=surrounding"

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
    --custom ${gnomad_del} \
    --no_stats \
    --fasta $FASTA_REF \
    --vcf \
    -o $OUTPUT_VCF_PATH

# generate tsv
echo -e "CHROM\tPOS\tID\tREF\tALT\t$(bcftools +split-vep -l $OUTPUT_VCF_PATH | cut -f2 | tr '\n' '\t' | sed 's/\t$//')" > $OUTPUT_TSV_PATH
bcftools +split-vep -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\n' -d -A tab $OUTPUT_VCF_PATH >> $OUTPUT_TSV_PATH

echo "$(date '+%Y-%m-%d %H:%M:%S') Job finished" >> ${logfile}
