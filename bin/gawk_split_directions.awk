## Split telomere file based on column 4 contents
## Date: 03/07/2025

BEGIN {
    FS="\t"; OFS="\t"
} {
    print > "direction."$3".telomere"
}
