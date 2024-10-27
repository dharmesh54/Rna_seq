#!/bin/bash

# Directory containing BAM files
BAM_DIR="/home/evpl/Dharmesh/rna_seq/bam_files"

# Loop through all BAM files in the directory
for bam_file in "$BAM_DIR"/*.bam; do
    # Define the output sorted BAM file name
    sorted_bam_file="${bam_file%.bam}_sorted.bam"
    
    # Sort the BAM file
    samtools sort -o "$sorted_bam_file" "$bam_file"
    
    # Print message indicating sorting completion
    echo "Sorted $bam_file to $sorted_bam_file"
done
