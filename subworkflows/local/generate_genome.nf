#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CUSTOM_GETCHROMSIZES  } from '../../modules/nf-core/custom/getchromsizes/main'
include { GNU_SORT              } from '../../modules/nf-core/gnu/sort'
include { GET_LARGEST_SCAFF     } from '../../modules/local/get_largest_scaff'

workflow GENERATE_GENOME {
    take:
    reference_file  // Channel: path(file)
    map_order       // Channel: val

    main:
    ch_versions     = Channel.empty()
    genome_size     = Channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //

    reference_file
        .combine(map_order)
        .map{ ref_meta, ref, map_order ->
             tuple(
                [   id : ref_meta,
                    map_order :map_order
                ],
                ref
             )
            }
        .branch{
            sorted      : it[0].map_order == "length"
            unsorted    : it[0].map_order == "unsorted"
        }
        .set{ch_genomesize_input}

    GENERATE_SORTED_GENOME (
        ch_genomesize_input.sorted
    )
    ch_versions         = ch_versions.mix( GENERATE_SORTED_GENOME.out.versions )
    ch_genomesize       = GENERATE_SORTED_GENOME.out.genomesize

    GENERATE_UNSORTED_GENOME (
        ch_genomesize_input.unsorted
    )
    ch_versions         = ch_versions.mix( GENERATE_UNSORTED_GENOME.out.versions )
    ch_genomesize       = ch_genomesize.mix(GENERATE_UNSORTED_GENOME.out.genomesize)

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
    dot_genome      = ch_genomesize
    ref_index       = CUSTOM_GETCHROMSIZES.out.fai
    versions        = ch_versions.ifEmpty(null)
}
