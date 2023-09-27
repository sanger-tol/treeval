#!/bin/bash
cram_path=$1
chunkn=0
for cram in ${cram_path}/*.cram; do
    rgline=$(samtools view -H $cram|grep "RG"|sed 's/\t/\\t/g'|sed "s/'//g")

    crampath=$(readlink -f ${cram})

    ncontainers=$(zcat ${crampath}.crai|wc -l)
    base=$(basename $cram .cram)

    from=0
    to=10000


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
done
