#!/usr/bin/env nextflow
//
// Check for synteny by aligning to fasta to reference genomes.
//
include { TEMP_MINIMAP2_ALIGN   } from '../../modules/local/TEMP_minimap_align'
include { GET_SYNTENY_GENOMES   } from '../../modules/local/get_synteny_genomes'

workflow SYNTENY {
    take:
        reference_tuple
        synteny_path
        assembly_classT

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: SEARCHES PREDETERMINED PATH FOR SYNTENIC GENOME FILES BASED ON CLASS
    //         EMITS PATH LIST
    //
    GET_SYNTENY_GENOMES(synteny_path, assembly_classT)

    //
    // LOGIC: GENERATES LIST OF GENOMES IN PATH AND BRANCHES ON WHETHER THERE IS DATA
    //
    GET_SYNTENY_GENOMES.out.genome_path
        .flatten()
        .branch { data -> 
            run: !data.toString().contains("empty")
            skip: data.toString().contains("empty")
        }
        .set { mm_intermediary }

    //
    // LOGIC: COMBINE WITH ABOVE .RUN CHANNEL ADD BOOLEANS FOR MINIMAP
    //
    reference_tuple
        .combine(mm_intermediary.run)
        .map { meta, fa, ref ->
            tuple([ id: meta.id,
                    single_end: true],
                fa, ref, false, false, true, false)
            }
        .set { mm_input }

    //
    // MODULE: ALIGNS THE SUNTENIC GENOMES TO THE REFERENCE GENOME
    //         EMITS ALIGNED PAF FILE
    //
    TEMP_MINIMAP2_ALIGN( mm_input.map { [it[0], it[1]] },
                    mm_input.map { it[2] },
                    mm_input.map { it[6] },
                    mm_input.map { it[3] },
                    mm_input.map { it[4] },
                    mm_input.map { it[5] }
    )
    ch_versions = TEMP_MINIMAP2_ALIGN.out.versions
    
    emit:
    ch_paf          = TEMP_MINIMAP2_ALIGN.out.paf
    versions        = ch_versions.ifEmpty(null)
}