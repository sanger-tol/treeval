#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { WINDOWMASKER_USTAT    } from '../../../modules/nf-core/windowmasker/ustat/main'
include { WINDOWMASKER_MKCOUNTS } from '../../../modules/nf-core/windowmasker/mkcounts/main'

workflow REPEAT_MASKING {
    take:
    ch_reference     // Channel: tuple [ val(meta), path(file) ]

    main:
    //
    // MODULE: MARK UP THE REPEAT REGIONS OF THE REFERENCE GENOME
    //
    WINDOWMASKER_MKCOUNTS (
        ch_reference
    )


    //
    // MODULE: CALCULATE THE STATISTICS OF THE MARKED UP REGIONS
    //
    WINDOWMASKER_USTAT(
        WINDOWMASKER_MKCOUNTS.out.counts,
        ch_reference
    )

    emit:
    repeat_intervals = WINDOWMASKER_USTAT.out.intervals
    _dummy           = channel.empty()  // To force the workflow to use named outputs under "out"
}
