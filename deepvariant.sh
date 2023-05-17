#!/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J DeepVariant      # Job name
#SBATCH -p ngs1gpu           # Partition Name
#SBATCH -c 6               # core preserved
#SBATCH --mem=90G           # memory used
#SBATCH --gres=gpu:1        # 使用的GPU數 請參考Queue資源設定
#SBATCH --mail-user=judychou60@gmail.com
#SBATCH --mail-type=END

wkdir=WKDIR
ID=SAMPLE_ID
ref_ver=hg38
file=${wkdir}/${ID}/${ID}_${ref_ver}_bwamem.markdup.recal.sorted.bam

# Reference path
HG38=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta

# Setup
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run_deepvariant_hg38.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

# environment setting
module load libs/singularity/3.7.1
workdir=${wkdir}/${ID}
cp ${HG38} -t ${workdir}
cp ${HG38}.fai -t ${workdir}

singularity run --nv -B ${workdir}:${workdir} \
    /opt/ohpc/Taiwania3/pkg/biology/DeepVariant/deepvariant_1.4.0-gpu.sif \
    /opt/deepvariant/bin/run_deepvariant \
    --ref=${workdir}/Homo_sapiens_assembly38.fasta \
    --model_type=WES \
    --reads=${file} \
    --output_vcf=${workdir}/${ID}.${ref_ver}.DeepVariant.vcf.gz \
    --output_gvcf=${workdir}/${ID}.${ref_ver}.DeepVariant.gvcf.gz \
    --num_shards=4