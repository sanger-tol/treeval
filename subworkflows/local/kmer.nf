#!/usr/bin/env nextflow

//
// Adapted from https://github.com/sanger-tol/genomeassembly
// the Sanger genomeassembly pipeline by @ksenia-krasheninnikova
//
// Convert BAM to CRAM, create index and calculate statistics
//

//
// MODULE IMPORT BLOCK
//
include { CAT_CAT             } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK         } from "../../modules/nf-core/fastk/fastk/main"
include { MERQURYFK_MERQURYFK } from '../../modules/nf-core/merquryfk/merquryfk/main'

workflow KMER {
    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    reads_path          // Channel: [ val(meta), val( str ) ]

    main:
    ch_versions                 = Channel.empty()

    //
    // LOGIC: PREPARE GET_READS_FROM_DIRECTORY INPUT
    //
    reference_tuple
        .combine( reads_path )
        .map { meta, ref, reads_path ->
            tuple(
                [   id          : meta.id,
                    single_end  : true  ],
                reads_path
            )
        }
        .set { get_reads_input }

    get_reads_input.view()

    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_read_paths       = GrabFiles( get_reads_input )

    ch_grabbed_read_paths.view()

    //
    // LOGIC: PACBIO READS FILES TO CHANNEL
    //
    // ch_grabbed_read_paths
    //     .map { meta, files ->
    //         tuple( files )
    //     }
    //     .flatten()
    //     .set { ch_reads }

    // ch_reads.view()
    // //
    // // LOGIC: 
    // //
    // reads_path
    //     .flatMap { meta, reads -> 
    //         reads instanceof List ? reads.collect{ [ meta, it ] } : [ [ meta, reads ] ] 
    //     }
    //     .set{ reads_ch }

    //
    // MODULE: 
    //
    CAT_CAT( ch_grabbed_read_paths )
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())

    //
    // LOGIC: 
    //
    CAT_CAT.out.file_out
        .map{ meta, reads -> 
            reads.getName().endsWith('gz') ? [meta, reads.getParent().toString() + '/' + reads.getBaseName().toString() + '.fa.gz'] : [meta, reads.getParent().toString() + '/' + reads.getBaseName().toString() + '.fa'] 
            }
        .set{ ch_reads_merged }

    //
    // LOGIC: 
    //
    CAT_CAT.out.file_out
        .join(ch_reads_merged)
        .map{ meta, reads_old, reads_new -> 
            reads_old.renameTo(reads_new); 
        }
    
    //
    // MODULE: 
    //
    FASTK_FASTK( ch_reads_merged )
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions.first())

    //
    // LOGIC: 
    //
    FASTK_FASTK.out.hist
        .join(FASTK_FASTK.out.ktab)
        .join(reference_tuple)
        .map{ meta, hist, ktab, meta_ref, primary -> 
                // hap.size() ? [ meta, hist, ktab, primary, hap ] :
                                [ meta, hist, ktab, primary, [] ] 
        } 
        .set{ ch_merq }

    ch_merq.view()
    //
    // MODULE: 
    //
    MERQURYFK_MERQURYFK ( ch_merq )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    emit:
    merquryk_completeness     = MERQURYFK_MERQURYFK.out.stats  // meta, stats
    merquryk_qv               = MERQURYFK_MERQURYFK.out.qv     // meta, qv
    versions                  = ch_versions.ifEmpty(null)
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.fasta.gz")

    "true"
}