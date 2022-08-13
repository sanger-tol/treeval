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
include { CAT_BLAST             } from '../../modules/local/cat_blast'
include { FILTER_BLAST          } from '../../modules/local/filter_blast'
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { PULL_DOT_AS           } from '../../modules/local/pull_dot_as'
include { GENERATE_GENOME       } from '../../modules/local/genome_generator'

workflow GENE_ALIGNMENT {
    ch_data             = Channel.value(params.alignment.geneset.toString())
                        .splitCsv()
                        .flatten()

    ch_datadir          = Channel.value(params.alignment.data_dir + params.assembly.class + '/csv_data/')

    SAMTOOLS_FAIDX ( [[params.assembly.sample], params.reference] )

    GENERATE_GENOME ( SAMTOOLS_FAIDX.out.fai )

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

    BLAST_BLASTN ( ch_alignment_data.others, BLAST_MAKEBLASTDB.out.db )

    // WAITING ON DECISION ON INCLUDING PROTEIN BLAST (BLASTX)

    grouped_tuple       = BLAST_BLASTN.out.txt
                        .map { meta, file ->
                            tuple([id: meta.org, type: meta.type],file) }
                        .groupTuple( by: [0] )

    CAT_BLAST ( grouped_tuple )

    FILTER_BLAST (CAT_BLAST.out.concat_blast)

    PULL_DOT_AS ( FILTER_BLAST.out.final_tsv )

}
