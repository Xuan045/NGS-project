#!/bin/bash

# Set the path to your CSV file
csv_file="/staging/biology/u4432941/SRA/sra_partae"

# Download settings
download_dir="/staging/biology/u4432941/SRA/download_files"
export PATH=$PATH:/work/opt/ohpc/Taiwania3/pkg/biology/SRAToolkit/sratoolkit_v3.0.3/bin

cd $download_dir
# Read the lines from the list
while IFS= read -r line; do
    # Extract the last column from the line
    sample_name=$(echo "$line" | awk -F',' '{print $NF}' | tr -d '\r\n')

    # Download files using SRAToolkit fastq-dump (download as fastq format and split into R1 and R2)
    fastq-dump --gzip --split-files "$sample_name"
    
    wait
done < $csv_file
