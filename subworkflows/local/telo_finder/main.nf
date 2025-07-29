#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GAWK as GAWK_UPPER_SEQUENCE   } from '../../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_REGIONS         } from '../../../modules/local/find/telomere_regions/main'
include { GAWK as GAWK_CLEAN_TELOMERE   } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_MAP_TELO         } from '../../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_WINDOWS         } from '../../../modules/local/find/telomere_windows/main'
include { EXTRACT_TELO                  } from '../../../modules/local/extract/telo/main'
include { TABIX_BGZIPTABIX              } from '../../../modules/nf-core/tabix/bgziptabix'

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
        FIND_TELOMERE_REGIONS.out.telomere
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_WINDOWS.out.versions )


    def windows_file = FIND_TELOMERE_WINDOWS.out.windows
    def fallback_file = GAWK_CLEAN_TELOMERE.out.output

    // Use EXTRACT_TELO if windows_file has content, otherwise fallback to GAWK_MAP_TELO
    def safe_windows = windows_file.ifEmpty { Channel.empty() }

    EXTRACT_TELO(
        safe_windows
    )
    ch_versions     = ch_versions.mix( EXTRACT_TELO.out.versions )

    GAWK_MAP_TELO(
        fallback_file,
        [],
        false
    )
    ch_versions     = ch_versions.mix( GAWK_MAP_TELO.out.versions )

    // Merge bed files into one for TABIX_BGZIPTABIX
    def merged_bed = EXTRACT_TELO.out.bed.mix(GAWK_MAP_TELO.out.output)

    TABIX_BGZIPTABIX(
        merged_bed
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )


    emit:
    bed_file        = EXTRACT_TELO.out.bed.ifEmpty { GAWK_MAP_TELO.out.output }
    bed_gz_tbi      = TABIX_BGZIPTABIX.out.gz_tbi
    bedgraph_file   = EXTRACT_TELO.out.bedgraph
    versions        = ch_versions
}
