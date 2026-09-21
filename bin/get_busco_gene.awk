## Get Busco genes
## Originally written by Yumi Sims (yy5)
## Re-written to use the GAWK module by Damon-Lee B Pointon (dp24)

BEGIN { OFS = "\t" }
$0 !~ /^#/ && $2 != "Missing" {

    if ($9 == "") {
        $9="no_orthodb_link"
    }

    sub(/:.*/, "", $3)

    start = ($4 < $5) ? $4 : $5
    end   = ($4 < $5) ? $5 : $4

    print $3, start, end, $1, $7, $6, $9
}
