#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { SEQTK_CUTN        } from '../../modules/nf-core/seqtk/cutn/main'
include { GAP_LENGTH        } from '../../modules/local/gap_length'
include { TABIX_BGZIPTABIX  } from '../../modules/nf-core/tabix/bgziptabix/main'

workflow GAP_FINDER {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(fasta) ]
    max_scaff_size      // Channel: val(size of largest scaffold in bp)

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
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAP_LENGTH (
        SEQTK_CUTN.out.bed
    )
    ch_versions     = ch_versions.mix( GAP_LENGTH.out.versions )

    //
    // LOGIC: Adding the largest scaffold size to the meta data so it can be used in the modules.config
    //
    SEQTK_CUTN.out.bed
        .combine(max_scaff_size)
        .map {meta, row, scaff ->
            tuple([ id          : meta.id,
                    max_scaff   : scaff
                ],
                file(row)
            )}
        .set { modified_bed_ch }

    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    TABIX_BGZIPTABIX (
        modified_bed_ch
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    gap_file        = GAP_LENGTH.out.bedgraph
    gap_tabix       = TABIX_BGZIPTABIX.out.gz_csi
    versions        = ch_versions.ifEmpty(null)
}
