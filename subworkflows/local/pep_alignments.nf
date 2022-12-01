#!/usr/bin/env nextflow

include { CONCAT_GFF            } from '../../modules/local/concat_gff'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/modules/nf-core/bedtools/sort/main'
include { TABIX_BGZIPTABIX      } from '../../modules/nf-core/modules/nf-core/tabix/bgziptabix/main'
include { MINIPROT_INDEX        } from '../../modules/sanger-tol/nf-core-modules/miniprot/index/main'
include { MINIPROT_ALIGN        } from '../../modules/sanger-tol/nf-core-modules/miniprot/align/main'

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
            index_file  : data[3]
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

    CONCAT_GFF ( grouped_tuple )

    BEDTOOLS_SORT ( CONCAT_GFF.out.concat_gff , 'gff')

    TABIX_BGZIPTABIX ( BEDTOOLS_SORT.out.sorted )

    emit:
    gff_file    = BEDTOOLS_SORT.out.sorted
    tbi_gff     = TABIX_BGZIPTABIX.out.gz_tbi
}