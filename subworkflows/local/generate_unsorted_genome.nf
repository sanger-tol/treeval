#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CUSTOM_GETCHROMSIZES  } from '../../modules/nf-core/custom/getchromsizes/main'

workflow GENERATE_UNSORTED_GENOME {
    take:
    reference_file  // Channel: path(file)

    main:
    ch_versions     = Channel.empty()
    genome_size     = Channel.empty()

    CUSTOM_GETCHROMSIZES (
        reference_file,
        "unsorted.genome"
        )
    ch_versions     = ch_versions.mix(  CUSTOM_GETCHROMSIZES.out.versions )

    emit:

    genomesize      = CUSTOM_GETCHROMSIZES.out.sizes
    ref_index       = CUSTOM_GETCHROMSIZES.out.fai
    versions        = ch_versions.ifEmpty(null)
}
