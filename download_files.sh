#!/bin/bash

# Set the path to your CSV file
csv_file="/staging/biology/u4432941/SRA/sra_result.csv"

# Read the CSV file line by line
while IFS= read -r line; do
    # Skip the header line
    if [[ $line == "Experiment Accession"* ]]; then
        continue
    fi

    # Extract the last column using awk
    sample_name=$(echo "$line" | awk -F',' '{print $NF}')

    # Do whatever you want with the sample_name variable
    echo "Sample name: $sample_name"

    # Add your additional logic here

done < "$csv_file"
