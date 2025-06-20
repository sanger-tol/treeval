## For lines !starting with '>' ensure they are upper characters
## Date: 17/06/2025

/^>/ {
    print
    next
}

{
    print toupper($0)
}
