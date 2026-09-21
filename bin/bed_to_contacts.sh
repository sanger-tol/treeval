#!/bin/bash

version='1.0.0'
if [ $1 == '-v' ];
then
    echo "$version"
else
    paste -d '\t' - - < $1 | awk 'BEGIN {FS="\t"; OFS="\t"} {if ($1 > $7) {print substr($4,1,length($4)-2),$12,$7,$8,"16",$6,$1,$2,"8",$11,$5} else {print substr($4,1,length($4)-2),$6,$1,$2,"8",$12,$7,$8,"16",$5,$11} }' | tr '\-+' '01'  | sort -k3,3d -k7,7d | awk 'NF==11'
fi
