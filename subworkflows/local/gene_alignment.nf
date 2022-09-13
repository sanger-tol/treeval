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
include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'

workflow GENE_ALIGNMENT {
    take:
    dot_genome // Channel: [val(meta), [ datafile ]]
    reference_tuple
    assembly_classT
    alignment_datadir
    alignment_genesets

    main:
    ch_versions         = Channel.empty()

    ch_data             = alignment_genesets
                            .splitCsv()
                            .flatten()

    ch_data
        .combine( alignment_datadir )
        .combine( assembly_classT )
        .set { csv_input } 

    CSV_GENERATOR ( csv_input.map { it[0] },
                    csv_input.map { it[1] },
                    csv_input.map { it[2] } )

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
    
    BLAST_MAKEBLASTDB ( reference_tuple.map { it[1] } )
    ch_versions         = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

    ch_alignment_data.others
        .combine( BLAST_MAKEBLASTDB.out.db )
        .multiMap { data ->
            input:          tuple (data[0],
                                    data[1]
                                )
            blast_db:       data[2]
        }
        .set{ blast_n_input }

    BLAST_BLASTN ( 
        blast_n_input.input,
        blast_n_input.blast_db
    )
    ch_versions         = ch_versions.mix(BLAST_BLASTN.out.versions)

    ch_alignment_data.pep
        .combine( BLAST_MAKEBLASTDB.out.db )
        .multiMap { data -> 
            input:          tuple (data[0],
                                    data[1]
                                )
            blast_db:       data[2]
        }
        .set{ blast_tn_input }

    BLAST_TBLASTN ( 
        blast_tn_input.input,
        blast_tn_input.blast_db
    )
    ch_versions         = ch_versions.mix(BLAST_TBLASTN.out.versions)

    BLAST_BLASTN.out.txt
        .mix(BLAST_TBLASTN.out.txt)
        .map { meta, file ->
            tuple([id: meta.org, type: meta.type], file) }
        .groupTuple( by: [0] )
        .set { grouped_tuple }

    CAT_BLAST ( grouped_tuple )

    FILTER_BLAST (CAT_BLAST.out.concat_blast) 
    ch_versions         = ch_versions.mix(FILTER_BLAST.out.versions)

    FILTER_BLAST.out.final_tsv
        .combine(dot_genome)
        .map { meta, tsv_file, org, genome ->
            tuple([ id          :   meta.id,
                    type        :   meta.type,
                    join_on     :   tsv_file.toString().split('-')[-2],
                    branch_by   :   tsv_file.toString().split('-')[-1].split('.tsv')[0]
            ],
            file(tsv_file), file(genome)
            )}
        .branch {
            meta, tsv, genome ->
            blast : meta.branch_by  == "filtered90"
                return [ meta, tsv, genome ]
            empty : meta.branch_by  == "EMPTY"
                return [ meta, tsv, genome ]
        }
        .set { bb_input }

    Channel
        .fromPath('../assets/gene_alignment/assm_*.as')
        .map { it -> 
            tuple ([ join_on    :   it.toString().split('/')[-1].split('_')[-1].split('.as')[0] ],
                    file(it)
                )}
        .set { as_file }

    bb_input.blast.join(as_file)

    bb_input.blast.view()

    UCSC_BEDTOBIGBED (
        bb_input.blast.map { [it[0], it[1]] },
        bb_input.blast.map { it[2] },
        bb_input.blast.map { it[3] })
    ch_versions         = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )
    
    emit:
    gene_alignment_bb   = UCSC_BEDTOBIGBED.out.bigbed
    versions            = ch_versions.ifEmpty(null)
}
