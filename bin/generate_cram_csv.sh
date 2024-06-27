#!/bin/bash

# generate_cram_csv.sh
# -------------------
# Generate a csv file describing the CRAM folder
# -------------------
# Author = yy5
# -------------------

cram_path=$1
chunkn=0
for cram in ${cram_path}/*.cram; do
    rgline=$(samtools view -H $cram | grep "@RG" -m1 | sed 's/\t/\\t/g' | sed -r 's/\\tDS.+$//')

    crampath=$(readlink -f ${cram})

    ncontainers=$(zcat ${crampath}.crai|wc -l)
    base=$(basename $cram .cram)

    from=0
    to=10000

    echo "START TO OUTPUT DATA FOR $crampath" 1>&2

    while [ $to -lt $ncontainers ]
    do
        echo $crampath,${crampath}.crai,${from},${to},${base},${chunkn},${rgline}
        from=$((to+1))
        ((to+=10000))
        ((chunkn++))
    done

    if [ $from -le $ncontainers ]
    then
        echo $crampath,${crampath}.crai,${from},${ncontainers},${base},${chunkn},${rgline}
        ((chunkn++))
    fi

    echo "CSV GEN COMPLETE FOR $crampath" 1>&2

done

echo "SCRIPT COMPLETE" 1>&2
