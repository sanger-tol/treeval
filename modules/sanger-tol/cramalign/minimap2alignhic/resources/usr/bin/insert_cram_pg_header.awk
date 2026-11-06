#!/usr/bin/awk -f
# Insert @PG lines from a saved CRAM header into a SAM/BAM stream.
#
# Usage: awk -v pgfile=<cram_pg.tmp> -f insert_cram_pg_header
#
# For every record stream (SAM/BAM piped via stdin), on the first @PG line:
#   1. Read saved CRAM @PG lines (pgfile) and heal provenance gaps:
#      if an entry lacks PP and immediately follows an entry that had PP,
#      link it to the preceding ID to close orphan branches (e.g. samtools -> lima.1).
#      Track all IDs to detect collisions.
#   2. If the current aligner's @PG ID collides with one in the CRAM history,
#      rename the aligner entry by appending .1 (.2, .3, …) until unique –
#      matching the samtools deduplication convention.
#   3. Strip any existing PP from the aligner line and inject a PP tag
#      linking it to the final ID of the healed CRAM chain.
# All subsequent lines are passed through unchanged.

function get_id(record) {
    n = split(record, t, "\t")
    for (i = 1; i <= n; i++)
        if (t[i] ~ /^ID:/)
            return substr(t[i], 4)
    return ""
}

function has_pp(record) {
    n_pp = split(record, t_pp, "\t")
    for (i_pp = 1; i_pp <= n_pp; i_pp++)
        if (t_pp[i_pp] ~ /^PP:/)
            return 1
    return 0
}

/^@PG/ && !pg_i {
    aid = get_id($0)
    pg_count = 0
    last_id = ""
    last_had_pp = 0
    while ((getline line < pgfile) > 0) {
        id = get_id(line)
        cur_has_pp = has_pp(line)
        if (!cur_has_pp && last_had_pp && last_id != "") {
            line = line "\tPP:" last_id
            cur_has_pp = 1
        }
        print line
        pg_lines[++pg_count] = line
        if (id != "") {
            seen[id] = 1
            last_id = id
            last_had_pp = cur_has_pp
        }
    }
    close(pgfile)
    pg_i = 1
    if (aid in seen) {
        sfx = 1
        while ((aid "." sfx) in seen) sfx++
        new_aid = aid "." sfx
        n = split($0, af, "\t"); out = ""
        for (j=1; j<=n; j++) { f = af[j]; if (f ~ /^ID:/) f = "ID:" new_aid; out = out (j > 1 ? "\t" : "") f }
        $0 = out
    }
    gsub(/\tPP:[^\t\r\n]+/, "", $0)
    if (last_id != "") { print $0 "\tPP:" last_id; next }
}
{ print }
