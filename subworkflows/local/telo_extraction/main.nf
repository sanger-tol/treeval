include { GAWK as GAWK_CLEAN_TELOMERE   } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_MAP_TELO         } from '../../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_WINDOWS         } from '../../../modules/local/find/telomere_windows/main'
include { EXTRACT_TELOMERE              } from '../../../modules/local/extract/telomere/main'
include { TABIX_BGZIPTABIX              } from '../../../modules/nf-core/tabix/bgziptabix'

workflow TELO_EXTRACTION {
    take:
    telomere_file   //tuple(meta, file)

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: GENERATES A WINDOWS FILE FROM THE ABOVE
    //
    FIND_TELOMERE_WINDOWS (
        telomere_file
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_WINDOWS.out.versions )


    def windows_file = FIND_TELOMERE_WINDOWS.out.windows
    def safe_windows = windows_file.ifEmpty { Channel.empty() }


    //
    // MODULE: Extract the telomere data from the FIND_TELOMERE
    //          file and reformat into bed
    //
    EXTRACT_TELOMERE(
        windows_file
    )
    ch_versions     = ch_versions.mix( EXTRACT_TELOMERE.out.versions )


    //
    // MODULE: REFORMAT THE WINDOWS FILE
    //
    GAWK_MAP_TELO(
        safe_windows,
        [],
        false
    )
    ch_gawk_output  = GAWK_MAP_TELO.out.output.ifEmpty( Channel.empty() )
    ch_versions     = ch_versions.mix( GAWK_MAP_TELO.out.versions )


    //
    // MODULE: Merge bed files into one for TABIX_BGZIPTABIX
    //
    // EXTRACT_TELO is the more important of the two, then we go to fallback, then just stop no point in running on empty file.
    def merged_bed = EXTRACT_TELOMERE.out.bed.ifEmpty { ch_gawk_output }

    TABIX_BGZIPTABIX(
        merged_bed
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )


    emit:
    bed_file        = merged_bed
    bed_gz_tbi      = TABIX_BGZIPTABIX.out.gz_tbi
    bedgraph_file   = EXTRACT_TELOMERE.out.bedgraph
    versions        = ch_versions

}
