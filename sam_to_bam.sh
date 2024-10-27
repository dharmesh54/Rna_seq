#!/bin/bash

# Directory containing SAM files
input_dir="/path/to/your/sam_files"  # replace with your actual directory path

# Loop through each SAM file in the directory
for sam_file in "$input_dir"/*.sam; do
    # Define output BAM file name by replacing the .sam extension with .bam
    bam_file="${sam_file%.sam}.bam"

    # Convert SAM to BAM using samtools
    echo "Converting $sam_file to $bam_file..."
    samtools view -Sb "$sam_file" > "$bam_file"

    # Check if the conversion was successful
    if [ $? -eq 0 ]; then
        echo "Successfully converted $sam_file to $bam_file"
    else
        echo "Error converting $sam_file"
    fi
done
