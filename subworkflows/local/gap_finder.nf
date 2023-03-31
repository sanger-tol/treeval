#!/usr/bin/env nextflow

include { SEQTK_CUTN    } from '../../modules/local/seqtk_cutn'
include { GAP_LENGTH    } from '../../modules/local/gap_length'

workflow GAP_FINDER {
    take:
    reference_tuple

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        reference_tuple
    )
    ch_versions = ch_versions.mix( SEQTK_CUTN.out.versions )


    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE
    //
    GAP_LENGTH (
        SEQTK_CUTN.out.bed
    )

    emit:
    gap_file        = GAP_LENGTH.out.bed
    versions        = ch_versions.ifEmpty(null)
}
