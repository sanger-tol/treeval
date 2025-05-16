#!/usr/bin/env nextflow

//
// Adapted from MicroFinder.v0.2
// by Tom Mathers
//

//
// MODULE IMPORT BLOCK
//

include { MINIPROT_ALIGN        } from '../../modules/nf-core/miniprot/align/main'
include { MINIPROT_INDEX        } from '../../modules/nf-core/miniprot/index/main'
include { SORT_FASTA            } from '../../modules/local/sort_fasta'
include { MICROFINDER_FILTER    } from '../../modules/local/microfinder_filter'

workflow MICROFINDER {
    take:
    reference_tuple     // Channel: [ val(meta), path(file) ]
    ouput_prefix          // Channel: [ val(meta), val( str ) ]
    scaffold_length_cutoff
    pep_files

    main:
    ch_versions             = Channel.empty()

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
    // MODULE: FILTER HITS
    //
    MICROFINDER_FILTER ( 
        MINIPROT_ALIGN.out.gff,
        reference_tuple,
        ouput_prefix,
        scaffold_length_cutoff
    )
    ch_versions         = ch_versions.mix( MICROFINDER_FILTER.out.versions )

    //
    // MODULE: REORDER ASSEMBLY
    //
    SORT_FASTA ( 
        MICROFINDER_FILTER.out.tsv,
        ouput_prefix,
        scaffold_length_cutoff
    )
    ch_versions     = ch_versions.mix(SORT_FASTA.out.versions)

    emit:
    SORT_FASTA.out.fa
    versions                = ch_versions.ifEmpty(null)
}

