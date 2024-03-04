#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GET_LARGEST_SCAFF             } from '../../modules/local/get_largest_scaff'
include { GENERATE_UNSORTED_GENOME      } from '../../subworkflows/local/generate_unsorted_genome'
include { GENERATE_SORTED_GENOME        } from '../../subworkflows/local/generate_sorted_genome'

workflow GENERATE_GENOME {
    take:
    reference_file  // Channel: path(file)
    map_order       // Channel: val

    main:
    ch_versions     = Channel.empty()
    ch_genomesize   = Channel.empty()
    ch_genome_fai   = Channel.empty()

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //

    reference_file
        .combine(map_order)
        .map{ ref_meta, ref, map_order ->
             tuple(
                [   id: ref_meta.id,
                    map_order :map_order
                ],
                ref
             )
            }
        .branch{
            sorted      : it[0].map_order == "length"
            unsorted    : it[0].map_order != "length"
        }
        .set{ch_genomesize_input}

    //
    // SUBWORKFLOW: GENERATE CHROMOSOME SIZES FILE RANKED BY LENGTH (DEFINED BY USER)
    //
    GENERATE_SORTED_GENOME (
        ch_genomesize_input.sorted
    )
    ch_versions         = ch_versions.mix( GENERATE_SORTED_GENOME.out.versions )
    ch_genomesize       = GENERATE_SORTED_GENOME.out.genomesize
    ch_genome_fai       = GENERATE_SORTED_GENOME.out.ref_index
    ch_versions         = GENERATE_SORTED_GENOME.out.versions

    //
    // SUBWORKFLOW: GENERATE UNSORTED CHROMOSOME SIZES FILE (DEFINED BY USER)
    //
    GENERATE_UNSORTED_GENOME (
        ch_genomesize_input.unsorted
    )
    ch_versions         = ch_versions.mix( GENERATE_UNSORTED_GENOME.out.versions )
    ch_genomesize       = ch_genomesize.mix( GENERATE_UNSORTED_GENOME.out.genomesize )
    ch_genome_fai       = ch_genome_fai.mix( GENERATE_UNSORTED_GENOME.out.ref_index )
    ch_versions         = GENERATE_UNSORTED_GENOME.out.versions

    //
    // MODULE: Cut out the largest scaffold size and use as comparator against 512MB
    //          This is the cut off for TABIX using tbi indexes
    //
    GET_LARGEST_SCAFF (
        ch_genomesize
    )
    ch_versions     = ch_versions.mix( GET_LARGEST_SCAFF.out.versions )

    emit:
    max_scaff_size  = GET_LARGEST_SCAFF.out.scaff_size.toInteger()
    dot_genome      = ch_genomesize
    ref_index       = ch_genome_fai
    versions        = ch_versions.ifEmpty(null)
}
