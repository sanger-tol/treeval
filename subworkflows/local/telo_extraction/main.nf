include { GAWK as GAWK_CLEAN_TELOMERE   } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_MAP_TELO         } from '../../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_WINDOWS         } from '../../../modules/local/find/telomere_windows/main'
include { EXTRACT_TELO                  } from '../../../modules/local/extract/telo/main'
include { TABIX_BGZIPTABIX              } from '../../../modules/nf-core/tabix/bgziptabix'

workflow TELO_EXTRACTION {
    take:
    telomere_file //tuple(meta, file)

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: CLEAN THE .TELOMERE FILE IF CONTAINS "you screwed up" ERROR MESSAGE
    //          (LIKELY WHEN USING LOWERCASE LETTERS OR BAD MOTIF)
    //          WORKS BE RETURNING LINES THAT START WITH '>'
    //
    GAWK_CLEAN_TELOMERE (
        telomere_file,
        [],
        false
    )
    ch_versions     = ch_versions.mix( GAWK_CLEAN_TELOMERE.out.versions )


    //
    // MODULE: GENERATES A WINDOWS FILE FROM THE ABOVE
    //
    FIND_TELOMERE_WINDOWS (
        telomere_file
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_WINDOWS.out.versions )


    def windows_file = FIND_TELOMERE_WINDOWS.out.windows
    def fallback_file = GAWK_CLEAN_TELOMERE.out.output

    // Use EXTRACT_TELO if windows_file has content, otherwise fallback to GAWK_MAP_TELO
    def safe_windows = windows_file.ifEmpty { Channel.empty() }
    def fallback_valid = fallback_file.ifEmpty { Channel.empty() }

    EXTRACT_TELO(
        safe_windows
    )
    ch_versions     = ch_versions.mix( EXTRACT_TELO.out.versions )

    GAWK_MAP_TELO(
        fallback_valid,
        [],
        false
    )
    ch_gawk_output  = GAWK_MAP_TELO.out.output.ifEmpty( Channel.empty() )
    ch_versions     = ch_versions.mix( GAWK_MAP_TELO.out.versions )

    //
    // MODULE: Merge bed files into one for TABIX_BGZIPTABIX
    //
    // EXTRACT_TELO is the more important of the two, then we go to fallback, then just stop no point in running on empty file.
    def merged_bed = EXTRACT_TELO.out.bed.ifEmpty { ch_gawk_output }
    //def merged_bed = EXTRACT_TELO.out.bed.mix(GAWK_MAP_TELO.out.output)


    TABIX_BGZIPTABIX(
        merged_bed
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    bed_file        = merged_bed
    bed_gz_tbi      = TABIX_BGZIPTABIX.out.gz_tbi
    bedgraph_file   = EXTRACT_TELO.out.bedgraph
    versions        = ch_versions

}
