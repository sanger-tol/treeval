#!/bin/bash

# Assign inputs
fasta="$1"
fai="$2"
newfai=""
# Check if the FASTA file ends with .fa
if [[ "$fasta" =~ \.fa$ ]]; then
    $newfai=$fai
# Check if the FASTA file ends with .fasta
elif [[ "$fasta" =~ \.fasta$ ]]; then
    # Create new FAI filename by removing the .fasta suffix
    mv "$fai" "${fasta}.fai"
else
    echo "File suffix is neither .fa nor .fasta."
fi

