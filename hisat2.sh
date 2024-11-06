#!/bin/bash

# Define the path to the HISAT2 genome index prefix
HISAT2_INDEX="/mnt/c/Users/dharm/Documents/Biostate/index_file/GRCm39_index" 

# Create an output directory for alignments if it doesn’t exist
mkdir -p sam_file

# Loop through all R1 gzipped files
for R1_FILE in *_R1*.fastq.gz; do
    # Define the corresponding R2 file
    R2_FILE="${R1_FILE/_R1/_R2}"

    # Extract sample name by removing file extensions and R1/R2 tag
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R[12].fastq.gz//')

    # Print the names of the files
    echo "Processing sample: $SAMPLE_NAME"
    echo "R1 file: $R1_FILE"
    echo "R2 file: $R2_FILE"

    # Run HISAT2 to align the reads
    echo "Running HISAT2 for $SAMPLE_NAME..."
    hisat2 -p 4 -x "$HISAT2_INDEX" \
        -1 "$R1_FILE" \
        -2 "$R2_FILE" \
        -S "sam_file/${SAMPLE_NAME}.sam"

    echo "$SAMPLE_NAME alignment complete."
done

echo "All samples aligned."


