#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'
include { CUSTOM_GETCHROMSIZES  } from '../../modules/nf-core/custom/getchromsizes/main'
include { GET_LARGEST_SCAFF     } from '../../modules/local/get_largest_scaff'

workflow GENERATE_GENOME {
    take:
    assembly_id     // Channel val(assembly_id)
    reference_file  // Channel [ val(meta), path(file) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: GENERATES A REFERENCE DATA TUPLE
    //
    reference_file
        .combine( assembly_id )
        .map { file, sample_id ->
            tuple ([id: sample_id],
                    file)
        }
        .set { to_chromsize }

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    CUSTOM_GETCHROMSIZES (
        to_chromsize
    )
    ch_versions     = ch_versions.mix(  CUSTOM_GETCHROMSIZES.out.versions )


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
    dot_genome      = CUSTOM_GETCHROMSIZES.out.sizes
    ref_index       = CUSTOM_GETCHROMSIZES.out.fai
    reference_tuple = to_chromsize
    versions        = ch_versions.ifEmpty(null)
}
