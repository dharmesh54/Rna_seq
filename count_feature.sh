#!/bin/bash

# Path to the annotation GTF file
GTF_FILE="/path/to/annotation_file.gtf"  # Update with the actual path

# Output directory for the counts
OUTPUT_DIR="./featurecounts_output"
mkdir -p "$OUTPUT_DIR"

# Loop through each BAM file in the current directory
for BAM_FILE in *.bam; do
    # Get the base name of the file (without the .bam extension)
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

    # Output file name for each sample
    OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_counts.txt"

    # Run featureCounts
    featureCounts -T 4 -p -t exon -g gene_id -a "$GTF_FILE" -o "$OUTPUT_FILE" "$BAM_FILE"

    echo "Counts generated for $SAMPLE_NAME in $OUTPUT_FILE"
done

echo "All BAM files have been processed."
