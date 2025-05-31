#!/bin/bash

# Load PLINK module
module load PLINK/1.90b6.17

# Base directory for your data
BASE_DIR="/carnegie/nobackup/scratch/tbellagio/gea_grene-net/lfmm_full"

# Loop through the range of indices
for i in {1..19}
do
    # Define input and output paths
    INPUT_FILE="${BASE_DIR}/clumping/input_clumping_bio${i}.txt"
    OUTPUT_PATH="${BASE_DIR}/clumping/output/output_clumping_bio${i}"

    # Run PLINK clumping
    plink \
    --bfile /carnegie/nobackup/scratch/tbellagio/gea_grene-net/leave_1_out/clumping/grenenet_og \
    --clump ${INPUT_FILE} \
    --clump-p1 1.5453657571674064e-08 \
    --clump-r2 0.2 \
    --clump-kb 250 \
    --out ${OUTPUT_PATH}

    # Check if PLINK ran successfully
    if [ $? -eq 0 ]; then
        echo "PLINK clumping successful for bio_${i}"
    else
        echo "Error in PLINK clumping for bio_${i}"
    fi
done
