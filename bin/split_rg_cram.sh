#!/bin/bash

input_cram=$1

# Extract the read groups
read_groups=$(samtools view -H $input_cram | grep '@RG' | awk '{for(i=1;i<=NF;i++){if($i ~ /^ID:/){print substr($i,4)}}}')

# Loop through each read group and create a separate CRAM file
for rg in $read_groups; do
    output_cram="output_${rg}.cram"
    samtools view -h -r $rg -o $output_cram $input_cram
    samtools index $output_cram
done
