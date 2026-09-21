#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SEQTK_CUTN                } from '../../../modules/nf-core/seqtk/cutn/main'
include { GAWK as GAWK_GAP_LENGTH   } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIPTABIX          } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow GAP_FINDER {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(fasta) ]

    main:

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        reference_tuple
    )

    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAWK_GAP_LENGTH (
        SEQTK_CUTN.out.bed,
        file("${projectDir}/bin/gawk_gap_length.awk"),
        false
    )

    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    TABIX_BGZIPTABIX (
        SEQTK_CUTN.out.bed
    )

    emit:
    gap_file        = GAWK_GAP_LENGTH.out.output
    gap_tabix       = TABIX_BGZIPTABIX.out.gz_index
}
