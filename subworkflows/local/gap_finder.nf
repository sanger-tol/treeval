#!/usr/bin/env nextflow

include { SEQTK_CUTN        } from '../../modules/nf-core/seqtk/cutn/main'
include { GAP_LENGTH        } from '../../modules/local/gap_length'
include { GET_LARGEST_SCAFF } from '../../modules/local/get_largest_scaff'
include { TABIX_BGZIPTABIX  } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow GAP_FINDER {
    take:
    reference_tuple     // Channel [ val(meta), path(fasta) ]
    dot_genome

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        reference_tuple
    )
    ch_versions     = ch_versions.mix( SEQTK_CUTN.out.versions )

    //
    // MODULE: Cut out the largest scaffold size and use as comparator against 512MB
    //          This is the cut off for TABIX using tbi indexes
    //
    GET_LARGEST_SCAFF ( dot_genome )

    // 
    // join bed from SEQTK_CUTN
    // with scaff_size from GET_LARGEST_SCAFF
    // map to create a new meta with scaff_size
    // branch based on scaff_size

    //
    // LOGIC: Adding the largest scaffold size to the meta data so it can be used in the modules.config
    //    
    SEQTK_CUTN.out.bed
        .combine(GET_LARGEST_SCAFF.out.scaff_size.toInteger())
        .map {meta, row, scaff -> 
            tuple([ id          : meta.id, 
                    max_scaff   : scaff >= 500000000 ? 'csi': ''
                ],
                file(row)
            )}
        .set { modified_bed_ch }

    modified_bed_ch.view()

    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAP_LENGTH (
        SEQTK_CUTN.out.bed
    )
    ch_versions     = ch_versions.mix( GAP_LENGTH.out.versions )

    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    TABIX_BGZIPTABIX (
        modified_bed_ch
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    gap_file        = GAP_LENGTH.out.bed
    gap_tabix       = TABIX_BGZIPTABIX.out.gz_csi
    versions        = ch_versions.ifEmpty(null)
}
