#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CUSTOM_GETCHROMSIZES  } from '../../modules/nf-core/custom/getchromsizes/main'
include { GNU_SORT              } from '../../modules/nf-core/gnu/sort'
include { GET_LARGEST_SCAFF     } from '../../modules/local/get_largest_scaff'

workflow GENERATE_GENOME {
    take:
    reference_file  // Channel path(file)

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    CUSTOM_GETCHROMSIZES (
        reference_file,
        "temp.genome"
    )
    ch_versions     = ch_versions.mix(  CUSTOM_GETCHROMSIZES.out.versions )

    //
    // MODULE: SORT CHROM SIZES BY CHOM SIZE NOT NAME
    //
    GNU_SORT (
        CUSTOM_GETCHROMSIZES.out.sizes
    )

    //
    // MODULE: Cut out the largest scaffold size and use as comparator against 512MB
    //          This is the cut off for TABIX using tbi indexes
    //
    GET_LARGEST_SCAFF (
        CUSTOM_GETCHROMSIZES.out.sizes
    )
    ch_versions     = ch_versions.mix( GET_LARGEST_SCAFF.out.versions )

    emit:
    max_scaff_size  = GET_LARGEST_SCAFF.out.scaff_size.toInteger()
    dot_genome      = GNU_SORT.out.sorted
    ref_index       = CUSTOM_GETCHROMSIZES.out.fai
    versions        = ch_versions.ifEmpty(null)
}
