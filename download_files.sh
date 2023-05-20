#!/bin/bash

# Set the path to your CSV file
csv_file="/staging/biology/u4432941/SRA/sra_result.csv"

# Download settings
download_dir="/staging/biology/u4432941/SRA/download_files"
export PATH=$PATH:/work/opt/ohpc/Taiwania3/pkg/biology/SRAToolkit/sratoolkit_v3.0.3/bin

# Read the CSV file and store each line in an array
mapfile -t lines < "$csv_file"

cd $download_dir
# Iterate over each line (excluding the header line)
for ((i = 1; i < ${#lines[@]}; i++)); do
    # Extract the last column from the line
    line=${lines[$i]}
    sample_name=$(echo "$line" | awk -F',' '{print $NF}' | tr -d '\r\n')

    # Download files using SRAToolkit fastq-dump (download as fastq format and split into R1 and R2)
    fastq-dump --split-files $sample_name
#    prefetch $sample_name
    wait

done

