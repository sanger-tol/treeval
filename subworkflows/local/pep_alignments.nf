#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CAT_CAT               } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/bedtools/sort/main'
include { TABIX_BGZIPTABIX      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { MINIPROT_INDEX        } from '../../modules/nf-core/miniprot/index/main'
include { MINIPROT_ALIGN        } from '../../modules/nf-core/miniprot/align/main'
include { GFF_TO_BED            } from '../../modules/local/gff_to_bed'

workflow PEP_ALIGNMENTS {
    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    pep_files           // Channel [ val(meta), path(file) ]
    max_scaff_size      // Channel val(size of largest scaffold in bp)

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
        .multiMap { data -> 
            pep_tuple   : tuple( [  id:     data[0].id,
                                    type:   data[0].type,
                                    org:    data[0].org
                                ],
                                data[1] )
            index_file  : tuple( [  id: "Reference",
                                ],
                                data[3] )
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
        .map { it ->
            tuple([ id:     it[0].org + '_pep',
                    type:   it[0].type
                ],
                it[1] )
        }
        .groupTuple( by: [0] )
        .set { grouped_tuple }

    //
    // MODULE: AS ABOVE OUTPUT IS BED FORMAT, IT IS MERGED PER ORGANISM + TYPE
    //
    CAT_CAT ( grouped_tuple )
    ch_versions         = ch_versions.mix( CAT_CAT.out.versions )

    //
    // MODULE: SORTS ABOVE OUTPUT AND RETAINS GFF SUFFIX
    //         EMITS A MERGED GFF FILE
    //
    BEDTOOLS_SORT ( CAT_CAT.out.file_out , [] )
    ch_versions         = ch_versions.mix( BEDTOOLS_SORT.out.versions )

    //
    // MODULE: CUTS GFF INTO PUNCHLIST
    //
    GFF_TO_BED ( CAT_CAT.out.file_out )
    ch_versions         = ch_versions.mix( GFF_TO_BED.out.versions )

    BEDTOOLS_SORT.out.sorted
        .combine(max_scaff_size)
        .map {meta, row, scaff -> 
            tuple([ id          : meta.id, 
                    max_scaff   : scaff >= 500000000 ? 'csi': ''
                ],
                file(row)
            )}
        .set { modified_bed_ch }

    //
    // MODULE: COMPRESS AND INDEX MERGED.GFF
    //         EMITS A TBI FILE
    //
    TABIX_BGZIPTABIX ( modified_bed_ch )
    ch_versions         = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    emit:
    gff_file            = BEDTOOLS_SORT.out.sorted
    tbi_gff             = TABIX_BGZIPTABIX.out.gz_tbi
    pep_punch           = GFF_TO_BED.out.punchlist
    versions            = ch_versions.ifEmpty(null)
}
