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
include { BB_GENERATOR          } from '../../modules/local/bb_generator.nf'

workflow GENE_ALIGNMENT {
    take:
    dot_genome // Channel: [val(meta), [ datafile ]]

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

    // WAITING ON DECISION ON INCLUDING PROTEIN BLAST (BLASTX)

    grouped_tuple       = BLAST_BLASTN.out.txt
                        .map { meta, file ->
                            tuple([id: meta.org, type: meta.type],file) }
                        .groupTuple( by: [0] )

    CAT_BLAST ( grouped_tuple )

    FILTER_BLAST (CAT_BLAST.out.concat_blast)
    ch_versions = ch_versions.mix(FILTER_BLAST.out.versions)

    PULL_DOT_AS ( FILTER_BLAST.out.final_tsv )

    blast               = FILTER_BLAST.out.final_tsv
    dotas               = PULL_DOT_AS.out.dotas
    blast_and_dotas     = blast.join(dotas)

    // Reformat blast_and_dotas, calculate value to branch on based on the TSV (returns filtered90 ?: EMPTY)
    blast_and_dotas
        .combine( dot_genome )
        .map { meta, tsv_file, as_file, source, genome_file ->
            tuple([ id          :   meta.id,
                    type        :   meta.type,
                    branch_by   :   tsv_file.toString().split('-')[-1].split('.tsv')[0]
            ],
            file(tsv_file), file(as_file), file(genome_file)
            )
        }
        .branch {
            meta, tsv, asm, geno ->
            blast : meta.branch_by  == "filtered90"
                return [ meta, tsv, asm, geno ]
            empty : meta.branch_by  == "EMPTY"
                return [ meta, tsv, asm, geno ]
        }
        .set { ch_todo }

    BB_GENERATOR ( ch_todo.blast )

    BB_GENERATOR.out.bb_out.view()

    emit:
    bb_files            = BB_GENERATOR.out.bb_out

    versions            = ch_versions.ifEmpty(null)
}
