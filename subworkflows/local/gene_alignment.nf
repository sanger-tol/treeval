#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { CSV_GENERATOR         } from '../../modules/local/csv_generator'
include { CSV_PULL              } from '../../modules/local/csv_pull'
include { PEP_ALIGNMENTS        } from './pep_alignments'
include { NUC_ALIGNMENTS        } from './nuc_alignments'

workflow GENE_ALIGNMENT {
    take:
    dot_genome          // Channel: [val(meta), [ datafile ]]
    reference_tuple
    assembly_classT
    alignment_datadir
    alignment_genesets
    alignment_common
    intron_size
    as_files

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: TAKES A SINGLE LIKE CSV STRING AND CONVERTS TO LIST OF VALUES
    //          LIST IS MERGED WITH DATA_DIRECTORY AND ORGANISM_CLASS
    //
    ch_data             = alignment_genesets
                            .splitCsv()
                            .flatten()

    ch_data
        .combine( alignment_datadir )
        .combine( assembly_classT )
        .set { csv_input } 
    //
    // MODULE: CONVERTS THE ABOVE VALUES INTO A STRING REPRESENTATIVE OF THE PATH
    //
    CSV_GENERATOR (     csv_input.map { it[0] },
                        csv_input.map { it[1] },
                        csv_input.map { it[2] } )

    //
    // MODULE: USES THE 
    //
    CSV_PULL (          csv_input.map { it[0] },
                        CSV_GENERATOR.out.csv_path)
    //
    // LOGIC: CONVERTS THE PATH OBJECT INTO A TUPLE OF
    //          [ [ META.ID, META.TYPE, META.ORG ], GENE_ALIGNMENT_FILE ]
    //          DATA IS THEN BRANCHED BASED ON META.TYPE TO THE APPROPRIATE
    //          SUBWORKFLOW
    //
    CSV_PULL.out.csv
        .splitCsv( header: true, sep:',')
        .map( row ->
        tuple([ org:    row.org,
                type:   row.type,
                id:     row.data_file.split('/')[-1].split('.MOD.')[0]
            ],
            file(row.data_file)
        ))
        .branch {
            pep: it[0].type  == 'pep'
            nuc: it[0].type  != 'pep'
        }
        .set {ch_alignment_data}

    pep_files = ch_alignment_data.pep.collect()
    nuc_files = ch_alignment_data.nuc.collect()

    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR PEPTIDE DATA, EMITS GFF AND TBI
    //
    PEP_ALIGNMENTS (    reference_tuple,
                        pep_files )
    
    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR RNA, NUCLEAR AND COMPLEMENT_DNA DATA, EMITS BIGBED
    //
    NUC_ALIGNMENTS (    reference_tuple,
                        nuc_files,
                        dot_genome,
                        intron_size )

    emit:
    pep_gff         = PEP_ALIGNMENTS.out.tbi_gff
    gff_file        = PEP_ALIGNMENTS.out.gff_file
    nuc_bb_files    = NUC_ALIGNMENTS.out.nuc_alignment
}
