#!/usr/bin/sh
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J IPAH_ANNOVAR         # Job name
#SBATCH -p ngs186G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 28               # 使用的core數 請參考Queue資源設定
#SBATCH --mem=186g           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o IPAH_ANNOVAR.out.log          # Path to the standard output file
#SBATCH -e IPAH_ANNOVAR.err.log          # Path to the standard error ouput file
#SBATCH --mail-user=    # email
#SBATCH --mail-type=END              # 指定送出email時機 可為NONE, BEGIN, END, FAIL, REQUEUE, ALL

### Please define the following variables
INPUT="/staging/biology/u4432941/SRA/output/ipah_bwamem.deepvariant.cohort.filtered.vcf.gz"
wkdir="/staging/biology/u4432941/SRA/output"
para="IPAH_deepvariant_ANNOVAR_hg38"

### DO NOT CHANGE
REF=/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg38/Homo_sapiens_assembly38.fasta
ANNOVAR="/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl"
humandb="/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/"

cd ${wkdir}
mkdir -p ${wkdir}
cd ${wkdir}
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${para}_run.log
exec >$logfile 2>&1
set -euo pipefail
set -x

bcftools norm -m- $INPUT -O z -o ${wkdir}/${para}.decom_hg38.vcf.gz
bcftools norm -f $REF ${wkdir}/${para}.decom_hg38.vcf.gz -O z -o ${wkdir}/${para}.decom_hg38.norm.vcf.gz

perl /opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/annotate_variation.pl ${wkdir}/${para}.decom_hg38.norm.vcf.gz -build hg38 -hgvs $humandb -out ${para}_hgvs
perl ${ANNOVAR} ${wkdir}/${para}.decom_hg38.norm.vcf.gz $humandb -buildver hg38 -out ${para} -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad30_genome,avsnp150,TWB1492_AF -operation gx,r,gx,gx,f,f,f  -arg '-hgvs',,,,,, -nastring . -vcfinput -polish

rm ${wkdir}/${para}.decom_hg38.vcf.gz ${para}.avinput
