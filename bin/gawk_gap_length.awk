## Calculate the length of a gap
## Date: 17/06/2025

BEGIN { OFS = "\t" }
{
    print $0, sqrt(($3-$2)*($3-$2))
}
