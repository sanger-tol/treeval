#!/bin/bash

# gff_to_bed.sh
# -------------------
# A shell script to convert a
# gff into bed format by stripping the
# PAF header and reorganising the data
# -------------------
# Author = yy5

version='1.0.0'

if [ $1 == '-v'];
then
    echo "$version"
else
    awk '{print $0"\t"sqrt(($3-$2)*($3-$2))}'|sed 's/\./0/g'
fi
