## Reformat the intercept file
## Date: 17/06/2025

function my_abs(x) {
    return x < 0 ? -x : x
}

{
    gsub(/\./, "0")
    printf "%s\t%.0f\n", $0, my_abs($3 - $2)
}
