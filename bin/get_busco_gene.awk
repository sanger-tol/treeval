## Get Busco genes
## Originally written by Yumi Sims (yy5)
## Re-written to use the GAWK module by Damon-Lee B Pointon (dp24)

BEGIN { OFS = "\t" }
$0 !~ /^#/ && $2 != "Missing" {

    if ($9 == "") {
        $9="no_orthodb_link"
    }

    sub(/:.*/, "", $3)

    if ($4 > $5) {
        print $3, $5, $4, $1, $7, $6, $9
    } else {
        print $3, $4, $5, $1, $7, $6, $9
    }
}
