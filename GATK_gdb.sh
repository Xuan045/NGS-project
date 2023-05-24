#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J NGS_GDB         # Job name
#SBATCH -p ngs53G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 8               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=53g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

wkdir=/staging/biology/u4432941/SRA/output
gdb_dir=${wkdir}/gatk_gdb
SampleList=${wkdir}/sample_list

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
logfile=./${TIME}_run_gatkGDB.log
exec > >(tee -a "$logfile") 2>&1
set -euo pipefail

##########################
# Consollidate GVCF
##########################
# Get sample name map
sample_name_map=${wkdir}/sample_name_map
while IFS= read -r line; do
    # Remove trailing carriage return (^M)
    ID=$(echo $line | tr -d '\r')

    cd ${wkdir}/${ID}
    id_path=$(realpath ${ID}_hg38_bwamem.HC.g.vcf.gz)
	echo -e "${ID}\t${id_path}" >> ${sample_name_map}

done<${SampleList}

${GATK4}/gatk GenomicsDBImport \
	--genomicsdb-workspace-path ${gdb_dir} \
	--sample-name-map ${sample_name_map}