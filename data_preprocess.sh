#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J NGS_GATK         # Job name
#SBATCH -p ngs92G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=92g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

wkdir=WKDIR
ID=SAMPLE_ID
R1=/staging/biology/u4432941/SRA/download_files/${ID}_1.fastq
R2=/staging/biology/u4432941/SRA/download_files/${ID}_2.fastq

# Reference path
HG38=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta
dbSNP=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/dbsnp_146.hg38.vcf.gz

# Tools
tool_dir=/opt/ohpc/Taiwania3/pkg/biology
GATK4=${tool_dir}/GATK/gatk_v4.2.3.0
PICARD=${tool_dir}/Picard/picard_v2.26.0/picard.jar
BWA=${tool_dir}/BWA/BWA_v0.7.17/bwa
SAMTOOLS=${tool_dir}/SAMTOOLS/samtools_v1.15.1/bin/samtools

# Setup
TIME=$(date +%Y%m%d%H%M)
logfile=./${TIME}_${ID}_run_preprocess.log
exec > >(tee -a "$logfile") 2>&1
set -euo pipefail
temp=$wkdir/$ID/temp
mkdir $temp

#####################
# BWA alignment
#####################
ref_ver=hg38
${BWA} mem -t 16 -R '@RG\tID:'${ID}'_'${ref_ver}'_bwamem\tLB:'${ID}'_'${ref_ver}'_bwamem\tSM:'${ID}'_'${ref_ver}'_bwamem\tPL:ILLUMINA\' ${HG38} ${R1} ${R2} > $temp/${ID}_${ref_ver}_bwamem.bam

java -Xmx92g -jar ${PICARD} SortSam \
    INPUT=$temp/${ID}_${ref_ver}_bwamem.bam \
    OUTPUT=$temp/${ID}_${ref_ver}_bwamem.sorted.bam \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

#####################
# Data preprocess
#####################
java -Xmx92g -jar ${PICARD} MarkDuplicates \
    INPUT=$temp/${ID}_${ref_ver}_bwamem.sorted.bam \
    OUTPUT=$temp/${ID}_${ref_ver}_bwamem.markdup.bam \
    METRICS_FILE=${ID}_${ref_ver}_bwamem.markdup.metrics \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

${GATK4}/gatk BaseRecalibrator \
    -R ${HG38} \
    --known-sites ${dbSNP} \
    -I $temp/${ID}_${ref_ver}_bwamem.markdup.bam \
    -O ${ID}_${ref_ver}_bwamem.markdup.recal_data.table \

${GATK4}/gatk ApplyBQSR \
    -R ${HG38} \
    -I $temp/${ID}_${ref_ver}_bwamem.markdup.bam \
    -bqsr ${ID}_${ref_ver}_bwamem.markdup.recal_data.table \
    -O $temp/${ID}_${ref_ver}_bwamem.markdup.recal.bam

java -Xmx92g -jar ${PICARD} SortSam \
    INPUT=$temp/${ID}_${ref_ver}_bwamem.markdup.recal.bam \
    OUTPUT=${ID}_${ref_ver}_bwamem.markdup.recal.sorted.bam \
    SORT_ORDER="coordinate" \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

##########################################
# Germline variant calling by GATK HC
##########################################
${GATK4}/gatk HaplotypeCaller \
    -R ${HG38} \
    -I ${ID}_${ref_ver}_bwamem.markdup.recal.sorted.bam \
    -O ${ID}_${ref_ver}_bwamem.HC.g.vcf.gz \
    -ERC GVCF

rm -rf $temp