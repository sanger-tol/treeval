#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CUSTOM_GETCHROMSIZES  } from '../../../modules/nf-core/custom/getchromsizes/main'
include { GNU_SORT              } from '../../../modules/nf-core/gnu/sort'

workflow GENERATE_SORTED_GENOME {
    take:
    reference_file  // Channel: path(file)

    main:
    ch_versions     = Channel.empty()
    genome_size     = Channel.empty()

    CUSTOM_GETCHROMSIZES (
        reference_file,
        "unsorted.genome"
        )
    ch_versions     = ch_versions.mix( CUSTOM_GETCHROMSIZES.out.versions )
    genome_size     = CUSTOM_GETCHROMSIZES.out.sizes

    GNU_SORT (
            CUSTOM_GETCHROMSIZES.out.sizes
            )
    ch_versions     = ch_versions.mix( GNU_SORT.out.versions )

    emit:
    genomesize      = GNU_SORT.out.sorted
    ref_index       = CUSTOM_GETCHROMSIZES.out.fai
    versions        = ch_versions
}
