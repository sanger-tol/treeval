#!/usr/bin/env nextflow

include { CAT_CAT               } from '../../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/bedtools/sort/main'
include { TABIX_BGZIPTABIX      } from '../../modules/nf-core/tabix/bgziptabix/main'
include { MINIPROT_INDEX        } from '../../modules/nf-core/miniprot/index/main'
include { MINIPROT_ALIGN        } from '../../modules/nf-core/miniprot/align/main'

workflow PEP_ALIGNMENTS {
    take:
    reference_tuple
    pep_files

    main:

    MINIPROT_INDEX ( reference_tuple )

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

    MINIPROT_ALIGN ( 
        formatted_input.pep_tuple,
        formatted_input.index_file
    )

    MINIPROT_ALIGN.out.gff
        .map { it ->
            tuple([ id:     it[0].org,
                    type:   it[0].type
                ],
                it[1] )
        }
        .groupTuple( by: [0] )
        .set { grouped_tuple }

    CAT_CAT ( grouped_tuple )

    BEDTOOLS_SORT ( CAT_CAT.out.file_out , 'gff')

    TABIX_BGZIPTABIX ( BEDTOOLS_SORT.out.sorted )

    emit:
    gff_file    = BEDTOOLS_SORT.out.sorted
    tbi_gff     = TABIX_BGZIPTABIX.out.gz_tbi
}
