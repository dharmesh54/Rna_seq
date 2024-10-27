#!/bin/bash

# Directory containing BAM files
input_dir="/path/to/your/bam_files"  # Replace with your actual directory path
threads=16  # Specify the number of threads to use

# Loop through each BAM file in the directory
for bam_file in "$input_dir"/*.bam; do
    # Define output BAM file name by appending "_no_duplicates" before the extension
    output_file="${bam_file%.bam}_no_duplicates.bam"
    metrics_file="${bam_file%.bam}_metrics.txt"

    # Run Picard to remove duplicates using multiple threads
    echo "Removing duplicates from $bam_file using $threads threads..."
    picard MarkDuplicates \
        I="$bam_file" \
        O="$output_file" \
        M="$metrics_file" \
        REMOVE_DUPLICATES=true \
        MAX_RECORDS_IN_RAM=5000000 \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORT_ORDER=coordinate \
        USE_JDK_DEFLATER=true \
        USE_JDK_INFLATER=true \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
        MAX_OPTICAL_DUPLICATE_SET_SIZE=300000 \
        READ_NAME_REGEX=null  # Adjust this based on your read name format if necessary

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully removed duplicates from $bam_file. Output saved as $output_file."
    else
        echo "Error removing duplicates from $bam_file."
    fi
done

