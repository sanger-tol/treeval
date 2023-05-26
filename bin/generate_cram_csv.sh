#!/bin/bash
cram=$1


chunkn=0

rgline=$(samtools view -H $cram|grep "RG"|sed 's/\t/\\t/g'|sed "s/'//g")

crampath=$(readlink -f ${cram})

ncontainers=$(zcat ${crampath}.crai|wc -l)
base=$(basename $cram .cram)

from=0
to=10000


while [ $to -le $ncontainers ]
do
    echo $crampath,${from},${to},${base},${chunkn},${rgline}
    from=$((to+1))
    ((to+=10000))
    ((chunkn++))
done

if [ $from -le $ncontainers ]
then
    
    echo $crampath,${from},${ncontainers},${base},${chunkn},${rgline}
    ((chunkn++))
fi
            
