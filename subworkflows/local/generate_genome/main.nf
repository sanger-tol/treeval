#!/usr/bin/env nextflow

//
// SUBWORKFLOW IMPORT BLOCK
//
include { GENERATE_UNSORTED_GENOME      } from '../generate_unsorted_genome/main'
include { GENERATE_SORTED_GENOME        } from '../generate_sorted_genome/main'

workflow GENERATE_GENOME {
    take:
    reference_file  // Channel: path(file)
    map_order       // Channel: val

    main:
    ch_genomesize   = channel.empty()
    ch_genome_fai   = channel.empty()


    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    reference_file
        .combine(map_order)
        .map{ ref_meta, ref, map_order_input ->
            tuple(
                [   id: ref_meta.id,
                    map_order :map_order_input
                ],
                ref
            )
        }
        .branch{ meta, _ref ->
            sorted      : meta.map_order == "length"
            unsorted    : meta.map_order != "length"
        }
        .set{ch_genomesize_input}


    //
    // SUBWORKFLOW: GENERATE CHROMOSOME SIZES FILE RANKED BY LENGTH (DEFINED BY USER)
    //
    GENERATE_SORTED_GENOME (
        ch_genomesize_input.sorted
    )
    ch_genomesize       = GENERATE_SORTED_GENOME.out.genomesize
    ch_genome_fai       = GENERATE_SORTED_GENOME.out.ref_index


    //
    // SUBWORKFLOW: GENERATE UNSORTED CHROMOSOME SIZES FILE (DEFINED BY USER)
    //
    GENERATE_UNSORTED_GENOME (
        ch_genomesize_input.unsorted
    )
    ch_genomesize       = ch_genomesize.mix( GENERATE_UNSORTED_GENOME.out.genomesize )
    ch_genome_fai       = ch_genome_fai.mix( GENERATE_UNSORTED_GENOME.out.ref_index )

    emit:
    dot_genome      = ch_genomesize
    ref_index       = ch_genome_fai
    ref             = reference_file
}
