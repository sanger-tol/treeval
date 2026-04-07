#!/usr/bin/env nextflow

//
// SUBWORKFLOW IMPORT BLOCK
//
include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx/main'
include { GNU_SORT               } from '../../../modules/nf-core/gnu/sort'

workflow GENERATE_GENOME {
    take:
    reference_file  // Channel: path(file)
    map_order       // Channel: val

    main:

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE MODIFIED FOR SCAFF SIZES
    //
    reference_file
        .combine(map_order)
        .map{ ref_meta, ref, map_order_input ->
            tuple(
                [   id: ref_meta.id,
                    map_order: map_order_input
                ],
                ref
            )
        }
        .set{ ch_genomesize_input }


    //
    // MODULE: INDEX THE INPUT FASTA
    //
    SAMTOOLS_FAIDX (
        ch_genomesize_input.map { meta, file -> [meta, file, []]},
        true // get sizes
    )


    //
    // MODULE: SORT THE INPUT FASTA FAI FOR LENGTH
    //
    GNU_SORT (
        SAMTOOLS_FAIDX.out.sizes.filter { meta, _ref -> meta.map_order == "length" }
    )


    //
    // LOGIC: IF length IS THE MAP ORDER, OUTPUT THE
    //        SORTED LENGTHS ELSE OUTPUT THE FAI
    //
    output = GNU_SORT.out.sorted.mix(
        SAMTOOLS_FAIDX.out.sizes.filter{ meta, _ref ->
            meta.map_order != "length"
        }
    )

    emit:
    dot_genome      = output
    ref_index       = SAMTOOLS_FAIDX.out.fai
    ref             = reference_file
}
