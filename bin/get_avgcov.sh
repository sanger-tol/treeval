#!/bin/bash

# get_avgcov.sh
# -------------------
# A shell script to calculate average coverage for each scaffold
# into bed format for use
# -------------------
# Author = yy5
# -------------------
version='1.0.0'
if [ $1 == '-v'];
then
    echo "$version"
else
    awk '{OFS="\t"; $5=$4*($3-$2); print}' $1|awk '{OFS="\t"; sum[$1]+=$5} END {for (chrom in sum) print chrom, sum[chrom]}'|awk 'BEGIN {FS="\t"; OFS="\t"} NR==FNR {genome[$1]=$2; next} {if ($1 in genome) print $1, genome[$1], $2, $3; else print $1, "NA", $2, $3}' -  $2| awk '{OFS="\t"; print $1,"0",$3,($2/$3)}' | awk 'BEGIN {FS="\t"; OFS="\t"} {printf "%s\t%s\t%s\t%.0f\n", $1, $2, $3, int($4 + 0.5)}'|sort -k1,1 -k2,2n> $3
fi