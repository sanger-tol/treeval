awk '{print $0"\t"sqrt(($3-$2)*($3-$2))}'|sed 's/\./0/g'
