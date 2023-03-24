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

    //
    // MODULE: CREATES INDEX OF REFERENCE FILE
    //
    MINIPROT_INDEX ( reference_tuple )

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

    //
    // MODULE: SORTS ABOVE OUTPUT AND RETAINS GFF SUFFIX
    //         EMITS A MERGED GFF FILE
    //
    BEDTOOLS_SORT ( CAT_CAT.out.file_out , [] )

    //
    // MODULE: COMPRESS AND INDEX MERGED.GFF
    //         EMITS A TBI FILE
    //
    TABIX_BGZIPTABIX ( BEDTOOLS_SORT.out.sorted )

    emit:
    gff_file    = BEDTOOLS_SORT.out.sorted
    tbi_gff     = TABIX_BGZIPTABIX.out.gz_tbi
}
