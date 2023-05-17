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

PIPELINE="/staging/biology/u4432941/SRA/data_preprocess.sh"
DAY=`date +%Y%m%d`

while read -r ID;
    do
    cd ${wkdir}
    mkdir -p ${ID}
    cd ${wkdir}/${ID}

    rsync ${PIPELINE} ./${DAY}_preprocess_${ID}.sh
    sed -i 's|WKDIR|'${wkdir}'|g' ./${DAY}_preprocess_${ID}.sh
    sed -i 's|SAMPLE_ID|'${ID}'|g' ./${DAY}_preprocess_${ID}.sh
    sbatch ./${DAY}_preprocess_${ID}.sh

    cd ${wkdir}

    done<${SampleList}
