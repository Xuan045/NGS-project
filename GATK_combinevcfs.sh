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
${GATK4}/gatk CombineGVCFs \
    -R ${HG38}
	--variant /staging/biology/u4432941/SRA/output/SRR18670381/SRR18670381_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670382/SRR18670382_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670383/SRR18670383_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670384/SRR18670384_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670385/SRR18670385_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670386/SRR18670386_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670387/SRR18670387_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670388/SRR18670388_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670389/SRR18670389_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670390/SRR18670390_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670391/SRR18670391_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670392/SRR18670392_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670393/SRR18670393_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670394/SRR18670394_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670395/SRR18670395_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670396/SRR18670396_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670397/SRR18670397_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670398/SRR18670398_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670399/SRR18670399_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670400/SRR18670400_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670401/SRR18670401_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670402/SRR18670402_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670403/SRR18670403_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670404/SRR18670404_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670405/SRR18670405_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670406/SRR18670406_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670407/SRR18670407_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670408/SRR18670408_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670409/SRR18670409_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670410/SRR18670410_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670411/SRR18670411_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670412/SRR18670412_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670413/SRR18670413_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670414/SRR18670414_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670415/SRR18670415_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670416/SRR18670416_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670417/SRR18670417_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670418/SRR18670418_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670419/SRR18670419_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670420/SRR18670420_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670421/SRR18670421_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670422/SRR18670422_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670423/SRR18670423_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670424/SRR18670424_hg38_bwamem.HC.g.vcf.gz \
    --variant /staging/biology/u4432941/SRA/output/SRR18670425/SRR18670425_hg38_bwamem.HC.g.vcf.gz \
	--O ${wkdir}/IPAH_cohort.g.vcf.gz