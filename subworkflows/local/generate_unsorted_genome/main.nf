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

    SAMTOOLS_FAIDX (
        ch_faidx_input,
        true // get sizes
    )

    SAMTOOLS_FAIDX.out.sizes
        .map { meta, sizes ->
            tuple( 
                meta, 
                file(sizes).renameTo("${meta.id}.unsorted.genome") 
                )
        }
        .set { ch_sizes }

    emit:

    genomesize      = ch_sizes
    ref_index       = SAMTOOLS_FAIDX.out.fai
    versions        = ch_versions
}
