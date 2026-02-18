#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx/main'

workflow GENERATE_UNSORTED_GENOME {
    take:
    reference_file  // Channel: path(file)

    main:
    ch_versions     = channel.empty()

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


    //
    // LOGIC: RENAME THE SORTED SIZES FILE
    //
    SAMTOOLS_FAIDX.out.sizes
        .map { meta, sizes ->
            tuple(
                meta,
                file(sizes).moveTo("${meta.id}.unsorted.genome")
                )
        }
        .set { ch_sizes }


    emit:

    genomesize      = ch_sizes
    ref_index       = SAMTOOLS_FAIDX.out.fai
    versions        = ch_versions
}
