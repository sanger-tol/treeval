#!/usr/bin/env nextflow

import java.math.RoundingMode;
import java.math.BigDecimal;

//
// MODULE IMPORT BLOCK
//
include { CAT_CAT               } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/bedtools/sort/main'
include { TABIX_BGZIPTABIX      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { MINIPROT_INDEX        } from '../../modules/nf-core/miniprot/index/main'
include { MINIPROT_ALIGN        } from '../../modules/nf-core/miniprot/align/main'
include { EXTRACT_COV_IDEN      } from '../../modules/local/extract_cov_iden'

workflow PEP_ALIGNMENTS {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    pep_files           // Channel: tuple [ val(meta), path(file) ]

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: CREATES INDEX OF REFERENCE FILE
    //
    MINIPROT_INDEX ( reference_tuple )
    ch_versions         = ch_versions.mix( MINIPROT_INDEX.out.versions )

    //
    // LOGIC: GETS LIST OF META AND PEP FILES FROM GENE_ALIGNMENT
    //        COMBINES WITH MINIPROT_INDEX OUTPUT
    //        CONVERTS TO TWO TUPLES FOR PEP DATA AND REFERENCE
    //
    pep_files
        .flatten()
        .buffer( size: 2 )
        .combine ( MINIPROT_INDEX.out.index )
        .multiMap { pep_meta, pep_file, miniprot_meta, miniprot_index ->
            pep_tuple   : tuple( [  id:     pep_meta.id,
                                    type:   pep_meta.type,
                                    org:    pep_meta.org
                                ],
                                pep_file )
            index_file  : tuple( [  id: "Reference",
                                ],
                                miniprot_index )
        }
        .set { formatted_input }

    //
    // MODULE: ALIGNS PEP DATA WITH REFERENCE INDEX
    //         EMITS GFF FILE
    //
    MINIPROT_ALIGN (
        formatted_input.pep_tuple,
        formatted_input.index_file
    )
    ch_versions         = ch_versions.mix( MINIPROT_ALIGN.out.versions )

    //
    // LOGIC: GROUPS OUTPUT GFFS BASED ON QUERY ORGANISMS AND DATA TYPE (PEP)
    //
    MINIPROT_ALIGN.out.gff
        .map { meta, file ->
            tuple(
                    [   id      :   meta.org + '_pep',
                        type    :   meta.type  ],
                    file
            )
        }
        .groupTuple( by: [0] )
        .set { grouped_tuple }

    //
    // MODULE: AS ABOVE OUTPUT IS BED FORMAT, IT IS MERGED PER ORGANISM + TYPE
    //
    CAT_CAT (
        grouped_tuple
    )
    ch_versions         = ch_versions.mix( CAT_CAT.out.versions )

    //
    // LOGIC: ADDING LINE COUNT TO THE FILE FOR BETTER RESOURCE USAGE
    //
    CAT_CAT.out.file_out
        .map { meta, file ->
            tuple ( [   id:     meta.id,
                        lines:  file.countLines()
                    ],
                    file
            )
        }
        .set { bedtools_input }

    //
    // MODULE: SORTS ABOVE OUTPUT AND RETAINS GFF SUFFIX
    //         EMITS A MERGED GFF FILE
    //
    BEDTOOLS_SORT (
        bedtools_input ,
        []
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_SORT.out.versions )

    //
    // MODULE: CUTS GFF INTO PUNCHLIST
    //
    EXTRACT_COV_IDEN (
        CAT_CAT.out.file_out
    )
    ch_versions         = ch_versions.mix( EXTRACT_COV_IDEN.out.versions )

    //
    // MODULE: COMPRESS AND INDEX MERGED.GFF
    //         EMITS A TBI FILE
    //
    TABIX_BGZIPTABIX (
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions         = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    gff_file            = BEDTOOLS_SORT.out.sorted
    tbi_gff             = TABIX_BGZIPTABIX.out.gz_tbi
    pep_punch           = EXTRACT_COV_IDEN.out.punchlist
    versions            = ch_versions.ifEmpty(null)
}
