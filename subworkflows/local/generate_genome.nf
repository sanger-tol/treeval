#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'
include { GENERATE_GENOME_FILE  } from '../../modules/local/generate_genome_file'
include { GET_LARGEST_SCAFF     } from '../../modules/local/get_largest_scaff'

workflow GENERATE_GENOME {
    take:
    assembly_id     // Channel val(assembly_id)
    reference_file  // Channel [ val(meta), path(file) ]

    main:
    ch_versions     = Channel.empty()

    //
    // LOGIC: GENERATES A REFERENCE DATA TUPLE
    //
    reference_file
        .combine( assembly_id )
        .map { it ->
            tuple ([id: it[1]],
                    it[0])
        }
        .set { to_samtools }

    //
    // MODULE: GENERATE INDEX OF REFERENCE
    //          EMITS REFERENCE INDEX FILE
    //
    SAMTOOLS_FAIDX ( to_samtools, [[],[]] )
    ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // MODULE: TRIMS INDEX INTO A GENOME DESCRIPTION FILE
    //         EMITS REFERENCE GEOME FILE AND REFERENCE INDEX FILE
    GENERATE_GENOME_FILE ( SAMTOOLS_FAIDX.out.fai )
    ch_versions     = ch_versions.mix( GENERATE_GENOME_FILE.out.versions )

    //
    // MODULE: Cut out the largest scaffold size and use as comparator against 512MB
    //          This is the cut off for TABIX using tbi indexes
    //
    GET_LARGEST_SCAFF ( GENERATE_GENOME_FILE.out.dotgenome )
    ch_versions     = ch_versions.mix( GET_LARGEST_SCAFF.out.versions )
 
    emit:
    max_scaff_size  = GET_LARGEST_SCAFF.out.scaff_size.toInteger()
    dot_genome      = GENERATE_GENOME_FILE.out.dotgenome
    ref_index       = SAMTOOLS_FAIDX.out.fai
    reference_tuple = to_samtools
    versions        = ch_versions.ifEmpty(null)
}
