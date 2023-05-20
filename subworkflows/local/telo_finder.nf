#!/usr/bin/env nextflow

include { FINDTELO                  } from '../../modules/local/findtelo'
include { FIND_TELOMERE_WINDOWS     } from '../../modules/local/find_telomere_windows'
include { EXTRACT_TELO              } from '../../modules/local/extract_telo'

workflow TELO_FINDER {

    take:
    reference_tuple     // Channel [ val(meta), path(fasta) ]
    teloseq

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATES TELOMERE TRACK FILES
    //
    FINDTELO (
        reference_tuple,
        teloseq
    )
    ch_versions     = ch_versions.mix( FINDTELO.out.versions )


    FIND_TELOMERE_WINDOWS (
        FINDTELO.out.telomere
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_WINDOWS.out.versions )

    EXTRACT_TELO (
        FIND_TELOMERE_WINDOWS.out.windows
    )
   ch_versions     = ch_versions.mix( EXTRACT_TELO.out.versions )

    emit:
    bed_file            = EXTRACT_TELO.out.bed
    bedgraph_file       = EXTRACT_TELO.out.bedgraph

    versions            = ch_versions.ifEmpty(null)
}
