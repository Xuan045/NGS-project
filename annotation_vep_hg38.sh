#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J VEP_hg38         # Job name
#SBATCH -p ngs186G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 28               # 使用的core數 請參考Queue資源設定 
#SBATCH --mem=186g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o $logfile          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=   # email
#SBATCH --mail-type=END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

source /work/opt/ohpc/Taiwania3/pkg/biology/Ensembl-VEP/Ensembl-VEP_v104.3/env.sh
export PATH=/opt/ohpc/Taiwania3/pkg/biology/HTSLIB/htslib_v1.13/bin/:$PATH
SOFTWARE=/work/opt/ohpc/Taiwania3/pkg/biology
REF=/staging/reserve/paylong_ntu/AI_SHARE/reference/VEP_ref
BCFTOOLS=/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools
VEP_CACHE_DIR=${SOFTWARE}/DATABASE/VEP/Cache
VEP_PLUGIN_DIR=${SOFTWARE}/DATABASE/VEP/Cache/Plugins
PLUGIN_DATA=${VEP_CACHE_DIR}/Data_for_plugins
VEP_PATH=${SOFTWARE}/Ensembl-VEP/Ensembl-VEP_v104.3
VEP_FASTA=${REF}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
ExACpLI_BAM=${REF}/GCF_000001405.38_GRCh38.p12_knownrefseq_alns.bam

INPUT_VCF_PATH=/staging/biology/u4432941/SRA/output/ipah_bwamem.deepvariant.cohort.filtered.vcf.gz
OUTPUT_VCF_PATH=/staging/biology/u4432941/SRA/output
SAMPLE_ID=IPAH_deepvariant_VEP_hg38

cd $OUTPUT_VCF_PATH
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${SAMPLE_ID}_run.log
exec >$logfile 2>&1
set -euo pipefail
set -x

$VEP_PATH/vep --cache --offline \
    --cache_version 104 \
    --merged \
    --assembly GRCh38 \
    --port 3337 \
    --dir_plugins $VEP_PLUGIN_DIR \
    --dir_cache $VEP_CACHE_DIR \
    -i $INPUT_VCF_PATH \
    --vcf \
    -o ${OUTPUT_VCF_PATH}/${SAMPLE_ID}.vcf \
    --check_existing \
    --plugin SpliceAI,snv=${PLUGIN_DATA}/spliceai_scores.raw.snv.hg38.vcf.gz,indel=${PLUGIN_DATA}/spliceai_scores.raw.indel.hg38.vcf.gz \
    --plugin LoFtool,${PLUGIN_DATA}/LoFtool_scores.txt \
    --plugin ExACpLI,${PLUGIN_DATA}/ExACpLI_values.txt \
    --fasta $VEP_FASTA \
    --bam $ExACpLI_BAM \
    --fork 20 \
    --force_overwrite \
    --use_given_ref

echo -e "CHROM\tPOS\tREF\tALT\t$($BCFTOOLS +split-vep -l ${OUTPUT_VCF_PATH}/${SAMPLE_ID}.vcf | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > ${wkdir}/IPAH_deepvariant_VEP_hg38.tsv
$BCFTOOLS +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\n' -d -A tab ${OUTPUT_VCF_PATH}/${SAMPLE_ID}.vcf >> ${wkdir}/IPAH_deepvariant_VEP_hg38.tsv
