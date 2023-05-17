#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J NGS_GATK         # Job name
#SBATCH -p ngs92G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 14               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=92g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=judychou60@gmail.com
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
GATK=${tool_dir}/GATK/gatk_v4.2.3.0
PICARD=${tool_dir}/Picard/picard_v2.26.0/picard.jar
BWA=${tool_dir}/BWA/BWA_v0.7.17/bwa
SAMTOOLS=${tool_dir}/SAMTOOLS/samtools_v1.15.1/bin/samtools

# Setup
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run_hg38.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail

#####################
# BWA alignment
#####################
ref_ver=hg38
${BWA} mem -t 16 -R '@RG\tID:'${ID}'_'${ref_ver}'_PL:ILLUMINA\' ${HG38} ${R1} ${R2} > ${ID}_${ref_ver}_bwamem.bam

java -Xmx48g -jar ${PICARD} SortSam \
    INPUT=${ID}_${ref_ver}_bwamem.bam \
    OUTPUT=${ID}_${ref_ver}_bwamem.sorted.bam \
    SORT_ORDER="coordinate" \
	CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

#####################
# Data preprocess
#####################
java -Xmx48g -jar ${PICARD} MarkDuplicates \
    INPUT=${ID}_${ref_ver}_bwamem.sorted.bam \
    OUTPUT=${ID}_${ref_ver}_bwamem.markdup.bam \
    METRICS_FILE=${ID}_${ref_ver}_bwamem.markdup.metrics \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

${GATK4}/gatk BaseRecalibrator \
    -R ${HG38} \
    -knownSites ${dbSNP} \
    -I ${ID}_${ref_ver}_bwamem.markdup.bam \
    -o ${ID}_${ref_ver}_bwamem.markdup.recal_data.grp \
    -rf BadCigar \
    -nct 16

${GATK4}/gatk PrintReads \
    -R ${HG38} \
    -I ${ID}_${ref_ver}_bwamem.markdup.bam \
    -BQSR ${ID}_${ref_ver}_bwamem.markdup.recal_data.grp \
    -o ${ID}_${ref_ver}_bwamem.markdup.recal.bam \
    -nct 16

java -Xmx48g -jar ${PICARD} SortSam \
    INPUT=${ID}_${ref_ver}_bwamem.markdup.recal.bam \
    OUTPUT=${ID}_${ref_ver}_bwamem.markdup.recal.sorted.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true

##########################################
# Germline variant calling by GATK HC
##########################################
${GATK4}/gatk HaplotypeCaller \
    -R ${HG38} \
    -I ${ID}_${ref_ver}_bwamem.markdup.recal.sorted.bam \
    -o ${ID}_${ref_ver}_bwamem.HC.g.vcf.gz \
    -ERC GVCF \
    --variant_index_type LINEAR \
    --variant_index_parameter 128000 \
    --dbsnp ${dbSNP} \
    --max_alternate_alleles 30 \
    -nct 16