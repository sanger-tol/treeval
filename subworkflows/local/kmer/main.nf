#!/usr/bin/env nextflow

//
// Adapted from https://github.com/sanger-tol/genomeassembly
// the Sanger genomeassembly pipeline by @ksenia-krasheninnikova
//
// Use FastK to count K-mers, plot spectra using MerquryFK
//

//
// MODULE IMPORT BLOCK
//
include { CAT_CAT             } from "../../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK         } from "../../../modules/nf-core/fastk/fastk/main"
include { MERQURYFK_MERQURYFK } from '../../../modules/nf-core/merquryfk/merquryfk/main'

workflow KMER {
    take:
    reference_tuple     // Channel: [ val(meta), path(file) ]
    reads_path          // Channel: [ val(meta), val( str ) ]

    main:

    //
    // LOGIC: PREPARE GET_READS_FROM_DIRECTORY INPUT
    //
    reads_path
        .map { meta, reads_path_input ->
            tuple(
                [   id          : meta.id,
                    single_end  : true  ],
                reads_path_input
            )
        }
        .set { get_reads_input }

    //
    // MODULE: JOIN PACBIO READ
    //
    CAT_CAT( get_reads_input )


    //
    // MODULE: COUNT KMERS
    //
    FASTK_FASTK( CAT_CAT.out.file_out )


    //
    // LOGIC: PREPARE MERQURYFK INPUT
    //
    FASTK_FASTK.out.hist
        .combine( FASTK_FASTK.out.ktab )
        .combine( reference_tuple )
        .map{ meta_hist, hist, _meta_ktab, ktab, _meta_ref, primary ->
            tuple( meta_hist, hist, ktab, primary, [] )
        }
        .set{ ch_merq }


    //
    // MODULE: USE KMER HISTOGRAM TO PRODUCE SPECTRA GRAPH
    //
    MERQURYFK_MERQURYFK (
        ch_merq,
        [[:], []],
        [[:], []]
    )

    emit:
    merquryk_completeness   = MERQURYFK_MERQURYFK.out.stats  // meta, stats
    merquryk_qv             = MERQURYFK_MERQURYFK.out.qv     // meta, qv
}
