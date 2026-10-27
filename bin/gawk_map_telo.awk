## Format the telo map file
## Date: 17/06/2025

BEGIN {
    OFS = "\t"
}

{
    sub(/^>/, "", $1)
    print $1, $4, $5, $6
}
