#!/bin/bash

# paf_to_bed12.sh
# -------------------
# A shell script to convert a
# paf into bed format for use
# in JBrowse
# -------------------
# Author = yy5

version='1.0.0'

if [ $1 == '-v'];
then
    echo "$version"
else
    cat $1 | awk 'BEGIN{FS="\t";}{a[$1]++;if(a[$1]==2)print v[$1] ORS $0;if(a[$1]>2)print;v[$1]=$0;}' | awk '$(NF+1) = ($10/$11)*100' | awk '$(NF+1) = ($10/$2)*100' | awk -vOFS='\t' '{print $6,$8,$9,$1,$2,$10,$(NF-1),$NF}' > $2
fi
