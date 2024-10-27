#!/bin/bash

#!/bin/bash

# Define the adapter file path (update this path to your adapter file)
ADAPTER_FILE="TruSeq3-PE-2.fa"  # Update with the correct adapter file path

# Loop through all R1 gzipped files
for R1_FILE in *_R1.fastq.gz; do
    # Define the corresponding R2 file based on R1 file name
    R2_FILE="${R1_FILE/_R1/_R2}"

    # Extract sample name by removing file extensions and R1/R2 tag
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1.fastq.gz//')

    # Print the names of the files
    echo "Processing sample: $SAMPLE_NAME"
    echo "R1 file: $R1_FILE"
    echo "R2 file: $R2_FILE"
    
    # Running Trimmomatic to remove adapters
    trimmomatic PE -threads 4 -phred33 \
        "$R1_FILE" \
        "$R2_FILE" \
        "${SAMPLE_NAME}_R1_paired.fq" \
        "${SAMPLE_NAME}_R1_unpaired.fq" \
        "${SAMPLE_NAME}_R2_paired.fq" \
        "${SAMPLE_NAME}_R2_unpaired.fq" \
        ILLUMINACLIP:"$ADAPTER_FILE":2:30:10

    echo "Trimming complete for $SAMPLE_NAME."
done

echo "All samples have been trimmed."



