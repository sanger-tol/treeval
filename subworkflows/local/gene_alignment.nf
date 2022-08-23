#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { CSV_GENERATOR         } from '../../modules/local/csv_generator'
include { BLAST_MAKEBLASTDB     } from '../../modules/nf-core/modules/blast/makeblastdb/main'
include { BLAST_BLASTN          } from '../../modules/nf-core/modules/blast/blastn/main'
include { BLAST_TBLASTN         } from '../../modules/sanger-tol/nf-core-modules/blast/tblastn/main'
include { CAT_BLAST             } from '../../modules/local/cat_blast'
include { FILTER_BLAST          } from '../../modules/local/filter_blast'
include { PULL_DOT_AS           } from '../../modules/local/pull_dot_as'
//include { BB_GENERATOR          } from '../../modules/local/bb_generator.nf'
include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'

workflow GENE_ALIGNMENT {
    take:
    dot_genome // Channel: [val(meta), [ datafile ]]

    main:
    ch_versions = Channel.empty()

    ch_data             = Channel.value(params.alignment.geneset.toString())
                        .splitCsv()
                        .flatten()

    ch_datadir          = Channel.value(params.alignment.data_dir + params.assembly.class + '/csv_data/')

    CSV_GENERATOR ( ch_datadir, ch_data )

    // Unique ID will be the org+chunk (size of the fasta for a dtype).
    CSV_GENERATOR.out.csv_path
        .splitCsv( header: true, sep:',')
        .map( row ->
        tuple([ org:    row.org,
                type:   row.type,
                id:     row.data_file.split('/')[-1].split('.MOD.')[0]
            ],
            file(row.data_file)
        ))
        .branch {
            pep: it[0].type     == 'pep'
            others: it[0].type  != 'pep'
            }
        .set {ch_alignment_data}

    BLAST_MAKEBLASTDB ( params.reference )
    ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    BLAST_BLASTN ( ch_alignment_data.others, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

    BLAST_TBLASTN ( ch_alignment_data.pep, BLAST_MAKEBLASTDB.out.db )
    ch_versions = ch_versions.mix(BLAST_TBLASTN.out.versions)

    BLAST_BLASTN.out.txt
        .mix(BLAST_TBLASTN.out.txt)
        .map { meta, file ->
            tuple([id: meta.org, type: meta.type], file) }
        .groupTuple( by: [0] )
        .set { grouped_tuple }

    //grouped_tuple       = BLAST_BLASTN.out.txt
    //                    .map { meta, file ->
    //                        tuple([id: meta.org, type: meta.type], file) }
    //                    .groupTuple( by: [0] )
    //

    grouped_tuple.view()

    CAT_BLAST ( grouped_tuple )

    FILTER_BLAST (CAT_BLAST.out.concat_blast) 
    ch_versions = ch_versions.mix(FILTER_BLAST.out.versions)

    blast               = FILTER_BLAST.out.final_tsv

    blast
        .map { meta, tsv_file ->
            tuple([ id          :   meta.id,
                    type        :   meta.type,
                    branch_by   :   tsv_file.toString().split('-')[-1].split('.tsv')[0]
            ],
            file(tsv_file)
            )}
        .branch {
            meta, tsv ->
            blast : meta.branch_by  == "filtered90"
                return [ meta, tsv ]
            empty : meta.branch_by  == "EMPTY"
                return [ meta, tsv ]
        }
        .set { ch_todo }

    ucsc_bb_input = ch_todo.blast
                        .combine(dot_genome)

    UCSC_BEDTOBIGBED (
          ucsc_bb_input.map { [it[0], it[1]] },
          ucsc_bb_input.map { it[3] }
    )
    
    emit:
    bb_files            = UCSC_BEDTOBIGBED.out.bigbed

    versions            = ch_versions.ifEmpty(null)
}
