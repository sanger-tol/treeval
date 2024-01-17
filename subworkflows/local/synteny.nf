#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'
include { GET_SYNTENY_GENOMES   } from '../../modules/local/get_synteny_genomes'

workflow SYNTENY {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    synteny_path        // Channel: val(meta)

    main:
    ch_versions                 = Channel.empty()

    reference_tuple
        .combine( synteny_path )
        .map { meta, reference, dir_path ->
            file("${dir_path}${meta.class}/*.fasta")
        }
        .set { syntenic_genomes }

    //
    // LOGIC: GENERATES LIST OF GENOMES IN PATH AND BRANCHES ON WHETHER THERE IS DATA
    //
    syntenic_genomes
        .flatten()
        .branch { data ->
            run                 : !data.toString().contains("empty")
            skip                : data.toString().contains("empty")
        }
        .set { mm_intermediary }

    //
    // LOGIC: COMBINE WITH ABOVE .RUN CHANNEL ADD BOOLEANS FOR MINIMAP
    //
    reference_tuple
        .combine( mm_intermediary.run )
        .multiMap { meta, syntenic_ref, ref ->
            syntenic_tuple  : tuple( meta, syntenic_ref )
            reference_fa    : ref
            bool_bam_output : false
            bool_cigar_paf  : true
            bool_cigar_bam  : false
        }
    .set { mm_input }

    //
    // MODULE: ALIGNS THE SUNTENIC GENOMES TO THE REFERENCE GENOME
    //         EMITS ALIGNED PAF FILE
    //
    MINIMAP2_ALIGN(
        mm_input.syntenic_tuple,
        mm_input.reference_fa,
        mm_input.bool_bam_output,
        mm_input.bool_cigar_paf,
        mm_input.bool_cigar_bam
    )
    ch_versions         = ch_versions.mix( MINIMAP2_ALIGN.out.versions )

    emit:
    ch_paf              = MINIMAP2_ALIGN.out.paf
    versions            = ch_versions.ifEmpty(null)
}
