#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J NGS_GATK         # Job name
#SBATCH -p ngs92G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=92g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=judychou60@gmail.com
#SBATCH --mail-type=END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

wkdir=/staging/biology/u4432941/SRA/output
merge_vcf=/staging/biology/u4432941/SRA/output/IPAH_cohort.g.vcf.gz

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
logfile=./${TIME}_run_gatkGC.log
exec > >(tee -a "$logfile") 2>&1
set -euo pipefail


###########################
# Perform joint genotyping
###########################
${GATK4}/gatk GenotypeGVCFs \
	-R ${HG38} \
	-V ${merge_vcf} \
	-O ${wkdir}/ipah_bwamem.gatkHC.merge.joint.vcf.gz

###################
# Filter variants
###################
${GATK4}/gatk VariantFiltration \
	-R ${HG38} \
	-V ${wkdir}/ipah_bwamem.gatkHC.merge.joint.vcf.gz \
	-O ${wkdir}/ipah_bwamem.gatkHC.merge.joint.filtered.vcf.gz \
	--filter-expression "DP < 5" \
	--filter-name "LowCoverage" \
	--filter-expression "DP >= 5 && DP < 30" \
	--filter-name "DepthLt30" \
	--filter-expression "QUAL < 30.0" \
	--filter-name "VeryLowQual" \
	--filter-expression "QUAL >= 30 && QUAL < 50.0" \
	--filter-name "LowQual" \
	--filter-expression "QD < 1.5" \
	--filter-name "LowQD"
