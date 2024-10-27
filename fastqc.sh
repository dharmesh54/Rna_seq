--threads#!/bin/bash

# Check if FastQC is installed
#if ! command -v fastqc &> /dev/null
#then
#echo "FastQC could not be found, please install FastQC."
#exit 1
#fi

# Loop through all .gz files in the current directory
for file in *.gz
do
    # Run FastQC on each file
    echo "Running FastQC on $file"
    fastqc --threads 2 "$file"
done

echo "FastQC analysis complete for all .gz files."
