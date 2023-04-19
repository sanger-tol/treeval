#!/bin/bash

# add_len2gap.sh
# ------------------------------------
# Script takes column two and three
# to calculate length of a gap in
# a given genomic fasta
# ------------------------------------
# Author = yy5

version="1.0.0"

if [ $1 == "-v" ]
then
    echo "$version"
else
    awk '{print $0"\t"sqrt(($3-$2)*($3-$2))}' $1
fi