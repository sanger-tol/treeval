#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'

workflow SYNTENY {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    synteny_path        // Channel: val(path)
    synteny_genomes     // Channel: val(name)

    main:
    ch_versions                 = Channel.empty()

    //
    // LOGIC: PULL SYNTENIC GENOMES FROM DIRECTORY STRUCTURE
    //          AND PARSE INTO CHANNEL PER GENOME
    //
    reference_tuple
        .combine( synteny_path )
        .combine( synteny_genomes )
        .map { meta, reference, dir_path, org ->
            file("${dir_path}${meta.class}/${org}.{fa,fasta}", checkIfExists: true)
        }
        .flatten()
        .combine( reference_tuple )
        .multiMap { syntenic_ref, meta, ref ->
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
