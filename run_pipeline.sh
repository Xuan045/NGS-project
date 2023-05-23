#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J RunPipeline         # Job name
#SBATCH -p ngs7G           # Partition Name 
#SBATCH -c 1               # core preserved 
#SBATCH --mem=7G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL

wkdir="/staging/biology/u4432941/SRA/output"
SampleList=${wkdir}/sample_list

pipeline=$1
PIPELINE="/staging/biology/u4432941/SRA/${pipeline}.sh"
DAY=`date +%Y%m%d`

while IFS= read -r line; do
    # Remove trailing carriage return (^M)
    ID=$(echo $line | tr -d '\r')

    cd ${wkdir}
    mkdir -p ${ID}
    cd ${wkdir}/${ID}

    rsync ${PIPELINE} ./${DAY}_${ID}_${pipeline}.sh
    sed -i 's|WKDIR|'${wkdir}'|g' ./${DAY}_${ID}_${pipeline}.sh
    sed -i 's|SAMPLE_ID|'${ID}'|g' ./${DAY}_${ID}_${pipeline}.sh
    sbatch ./${DAY}_${ID}_${pipeline}.sh

done<${SampleList}
