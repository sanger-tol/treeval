#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GAWK as GAWK_UPPER_SEQUENCE   } from '../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_REGIONS         } from '../../modules/local/find_telomere_regions'
include { GAWK as GAWK_CLEAN_TELOMERE   } from '../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_WINDOWS         } from '../../modules/local/find_telomere_windows'
include { EXTRACT_TELO                  } from '../../modules/local/extract_telo'
include { TABIX_BGZIPTABIX              } from '../../modules/nf-core/tabix/bgziptabix'

workflow TELO_FINDER {

    take:
    reference_tuple     // Channel: tuple [ val(meta), path(fasta) ]
    teloseq

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: UPPERCASE THE REFERENCE SEQUENCE
    //
    GAWK_UPPER_SEQUENCE(
        reference_tuple,
        [],
        false,
    )
    ch_versions     = ch_versions.mix( GAWK_UPPER_SEQUENCE.out.versions )

    //
    // MODULE: FINDS THE TELOMERIC SEQEUNCE IN REFERENCE
    //
    FIND_TELOMERE_REGIONS (
        GAWK_UPPER_SEQUENCE.out.output,
        teloseq
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_REGIONS.out.versions )


    //
    // MODULE: CLEAN THE .TELOMERE FILE IF CONTAINS "you screwed up" ERROR MESSAGE
    //          (LIKELY WHEN USING LOWERCASE LETTERS OR BAD MOTIF)
    //          WORKS BE RETURNING LINES THAT START WITH '>'
    //
    GAWK_CLEAN_TELOMERE (
        FIND_TELOMERE_REGIONS.out.telomere,
        [],
        false
    )
    ch_versions     = ch_versions.mix( GAWK_CLEAN_TELOMERE.out.versions )


    //
    // MODULE: GENERATES A WINDOWS FILE FROM THE ABOVE
    //
    FIND_TELOMERE_WINDOWS (
        GAWK_CLEAN_TELOMERE.out.output
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_WINDOWS.out.versions )

    //
    // MODULE: EXTRACTS THE LOCATION OF TELOMERIC SEQUENCE BASED ON THE WINDOWS
    //
    EXTRACT_TELO (
        FIND_TELOMERE_WINDOWS.out.windows
    )
    ch_versions     = ch_versions.mix( EXTRACT_TELO.out.versions )

    //
    // MODULE: BGZIP AND TABIX THE OUTPUT FILE
    //
    TABIX_BGZIPTABIX (
        EXTRACT_TELO.out.bed
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    bed_file        = EXTRACT_TELO.out.bed
    bed_gz_tbi      = TABIX_BGZIPTABIX.out.gz_tbi
    bedgraph_file   = EXTRACT_TELO.out.bedgraph
    versions        = ch_versions.ifEmpty(null)
}
