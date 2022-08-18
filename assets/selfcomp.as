table selfcomp
"selfcomp"
(
string  chrom;		"Reference sequence chromosome or scaffold id"
uint    chromStart;	"Start position of feature on chromosome"
uint    chromEnd;	"End position of feature on chromosome"
string  name;		"Querysequence chromosome or scaffold id"
string  length;     "Length of the motif"
char[1] strand;		"+ or - for strand of query seq"
string    qStart;     "Start position of feature on query sequence"
string    qEnd;       "End position of feature on query sequence"
string  qLength;    "Length of feature on query sequence"
)
