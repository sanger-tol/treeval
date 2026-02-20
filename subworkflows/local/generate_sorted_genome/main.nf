#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx/main'
include { GNU_SORT              } from '../../../modules/nf-core/gnu/sort'

workflow GENERATE_SORTED_GENOME {
    take:
    reference_file  // Channel: path(file)

    main:

    reference_file
        .map { ref_meta, ref ->
            tuple( ref_meta, ref, [] )
        }
        .set { ch_faidx_input }


    //
    // MODULE: INDEX THE INPUT FASTA
    //
    SAMTOOLS_FAIDX (
        ch_faidx_input,
        true // get sizes
    )

    SAMTOOLS_FAIDX.out.sizes
        .view{"SIZES: $it"}

    //
    // MODULE: SORT THE SIZES FILE
    //
    GNU_SORT (
        SAMTOOLS_FAIDX.out.sizes
    )


    //
    // LOGIC: RENAME THE SORTED SIZES FILE
    //
    // GNU_SORT.out.sorted
    //     .map { meta, sizes ->
    //         tuple(
    //             meta,
    //             file(sizes).copyTo("${meta.id}.sorted.genome")
    //             )
    //     }
    //     .set { ch_sizes }

    emit:
    genomesize      = GNU_SORT.out.sorted
    ref_index       = SAMTOOLS_FAIDX.out.fai
}
