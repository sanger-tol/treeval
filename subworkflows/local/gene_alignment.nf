#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { CSV_GENERATOR         } from '../../modules/local/csv_generator'

workflow GENE_ALIGNMENT {
    ch_data             = Channel.value(params.alignment.geneset.toString())
                        .splitCsv()
                        .flatten()

    ch_datadir          = Channel.value(params.alignment.data_dir + params.assembly.class + '/csv_data/')

    CSV_GENERATOR ( ch_datadir, ch_data )

}
