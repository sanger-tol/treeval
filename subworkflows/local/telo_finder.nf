#!/usr/bin/env nextflow

include { FINDTELO                  } from '../../modules/local/findtelo'
include { FIND_TELOMERE_WINDOWS     } from '../../modules/local/find_telomere_windows'

workflow TELO_FINDER {

    take:
    reference_tuple     // Channel [ val(meta), path(fasta) ]
    teloseq

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
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
   
    emit:
    windows_file        = FIND_TELOMERE_WINDOWS.out.windows
    versions            = ch_versions.ifEmpty(null)
}
