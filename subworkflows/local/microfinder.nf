#!/usr/bin/env nextflow

//
// Adapted from MicroFinder.v0.2
// by Tom Mathers
//
// Use FastK to count K-mers, plot spectra using MerquryFK
//

//
// MODULE IMPORT BLOCK
//

include { MINIPROT_ALIGN        } from '../../modules/nf-core/miniprot/align/main'

workflow MICROFINDER {
    take:
    reference_tuple     // Channel: [ val(meta), path(file) ]
    ouput_prefix          // Channel: [ val(meta), val( str ) ]
    scaffold_length_cutoff

    main:
    ch_versions             = Channel.empty()



    //
    // MODULE: ALIGNS PEP DATA WITH REFERENCE INDEX
    //         EMITS GFF FILE
    //
    MINIPROT_ALIGN (
        formatted_input.pep_tuple,
        formatted_input.index_file
    )
    ch_versions         = ch_versions.mix( MINIPROT_ALIGN.out.versions )


    // //
    // // MODULE: JOIN PACBIO READ
    // //
    // CAT_CAT( ch_grabbed_read_paths )
    // ch_versions             = ch_versions.mix( CAT_CAT.out.versions.first() )

    // //
    // // MODULE: COUNT KMERS
    // //
    // FASTK_FASTK( CAT_CAT.out.file_out )
    // ch_versions             = ch_versions.mix( FASTK_FASTK.out.versions.first() )

    // //
    // // LOGIC: PREPARE MERQURYFK INPUT
    // //
    // FASTK_FASTK.out.hist
    //     .combine( FASTK_FASTK.out.ktab )
    //     .combine( reference_tuple )
    //     .map{ meta_hist, hist, meta_ktab, ktab, meta_ref, primary ->
    //         tuple( meta_hist, hist, ktab, primary, [] )
    //     }
    //     .set{ ch_merq }

    // //
    // // MODULE: USE KMER HISTOGRAM TO PRODUCE SPECTRA GRAPH
    // //
    // MERQURYFK_MERQURYFK (
    //     ch_merq,
    //     [],
    //     []
    // )
    // ch_versions             = ch_versions.mix( MERQURYFK_MERQURYFK.out.versions.first() )

    emit:
    merquryk_completeness   = MERQURYFK_MERQURYFK.out.stats  // meta, stats
    merquryk_qv             = MERQURYFK_MERQURYFK.out.qv     // meta, qv
    versions                = ch_versions.ifEmpty(null)
}

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val( meta ), path( "in" )

    output:
    tuple val( meta ), path( "in/*.fasta.gz" )

    "true"
}
