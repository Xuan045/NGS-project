#!/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J CombineDeepVariant      # Job name
#SBATCH -p ngs1gpu           # Partition Name
#SBATCH -c 6               # core preserved
#SBATCH --mem=90G           # memory used
#SBATCH --gres=gpu:1        # 使用的GPU數 請參考Queue資源設定
#SBATCH --mail-user=judychou60@gmail.com
#SBATCH --mail-type=END

wkdir=/staging/biology/u4432941/SRA/output

# Reference path
HG38=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta

# Setup
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run_glnexus_hg38.log
exec > >(tee -a "$logfile") 2>&1
set -euo pipefail
set -x

# environment setting
module load libs/singularity/3.10.2
cp ${HG38} -t ${wkdir}
cp ${HG38}.fai -t ${wkdir}

singularity run --nv -B ${wkdir}:${wkdir} \
    /opt/ohpc/Taiwania3/pkg/biology/GLnexus/glnexus_v1.4.1.sif \
    /usr/local/bin/glnexus_cli \
    --config DeepVariantWGS \
    /staging/biology/u4432941/SRA/output/SRR18670381/SRR18670381.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670382/SRR18670382.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670383/SRR18670383.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670384/SRR18670384.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670385/SRR18670385.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670386/SRR18670386.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670387/SRR18670387.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670388/SRR18670388.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670389/SRR18670389.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670390/SRR18670390.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670391/SRR18670391.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670392/SRR18670392.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670393/SRR18670393.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670394/SRR18670394.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670395/SRR18670395.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670396/SRR18670396.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670397/SRR18670397.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670398/SRR18670398.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670399/SRR18670399.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670400/SRR18670400.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670401/SRR18670401.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670402/SRR18670402.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670403/SRR18670403.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670404/SRR18670404.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670405/SRR18670405.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670406/SRR18670406.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670407/SRR18670407.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670408/SRR18670408.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670409/SRR18670409.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670410/SRR18670410.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670411/SRR18670411.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670412/SRR18670412.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670413/SRR18670413.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670414/SRR18670414.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670415/SRR18670415.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670416/SRR18670416.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670417/SRR18670417.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670418/SRR18670418.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670419/SRR18670419.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670420/SRR18670420.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670421/SRR18670421.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670422/SRR18670422.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670423/SRR18670423.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670424/SRR18670424.hg38.DeepVariant.gvcf.gz \
    /staging/biology/u4432941/SRA/output/SRR18670425/SRR18670425.hg38.DeepVariant.gvcf.gz \
    | bcftools view - | bgzip -c > ${wkdir}/deepvariant.cohort.vcf.gz