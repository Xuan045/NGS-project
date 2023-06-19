#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J DeepVariant_filter         # Job name
#SBATCH -p ngs92G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=92g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

wkdir=/staging/biology/u4432941/SRA/output
vcf=/staging/biology/u4432941/SRA/output/ipah_bwamem.deepvariant.cohort.vcf.gz

# Reference path
HG38=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta

# Tools
tool_dir=/opt/ohpc/Taiwania3/pkg/biology
GATK4=${tool_dir}/GATK/gatk_v4.2.3.0
PICARD=${tool_dir}/Picard/picard_v2.26.0/picard.jar
BWA=${tool_dir}/BWA/BWA_v0.7.17/bwa
SAMTOOLS=${tool_dir}/SAMTOOLS/samtools_v1.15.1/bin/samtools

# Setup
TIME=$(date +%Y%m%d%H%M)
logfile=./${TIME}_deepvariant_filter.log
exec > >(tee -a "$logfile") 2>&1
set -euo pipefail

###################
# Filter variants
###################
bcftools norm -m- -Oz $vcf > ${wkdir}/ipah_bwamem.deepvariant.cohort.norm.vcf.gz
${GATK4}/gatk VariantFiltration \
	-R ${HG38} \
	-V ${wkdir}/ipah_bwamem.deepvariant.cohort.norm.vcf.gz $vcf \
	-O ${wkdir}/ipah_bwamem.deepvariant.cohort.norm.filtered.vcf.gz \
	--filter-expression "QUAL < 30" \
	--filter-name "LowQual" \
