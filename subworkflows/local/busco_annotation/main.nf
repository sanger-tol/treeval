#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

//
// MODULE IMPORT BLOCK
//
include { BUSCO_BUSCO                       } from '../../../modules/nf-core/busco/busco/main'
include { UCSC_BEDTOBIGBED                  } from '../../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                     } from '../../../modules/nf-core/bedtools/sort/main'
include { GAWK as GAWK_EXTRACT_BUSCOGENE    } from '../../../modules/nf-core/gawk/main'

//
// SUBWORKFLOW IMPORT BLOCK
//
include { ANCESTRAL_GENE                } from '../ancestral_gene/main'

workflow BUSCO_ANNOTATION {
    take:
    dot_genome           // Channel: tuple [val(meta), [ datafile ]]
    reference_tuple      // Channel: tuple [val(meta), [ datafile ]]
    lineageinfo          // Channel: val(lineage_db)
    lineagespath         // Channel: val(/path/to/buscoDB)
    buscogene_as         // Channel: val(dot_as location)
    ancestral_table      // Channel: val(ancestral_table location)

    main:
    ch_versions                 = Channel.empty()

    //
    // MODULE: RUN BUSCO TO OBTAIN FULL_TABLE.CSV
    //         EMITS FULL_TABLE.CSV
    //
    BUSCO_BUSCO (
        reference_tuple,
        "genome",
        lineageinfo,
        lineagespath,
        []
    )
    ch_versions          = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    ch_busco_full_table  = BUSCO_BUSCO.out.busco_dir.map { meta, dir -> tuple(meta, files(dir.resolve("*/*/full_table.tsv"), checkIfExists: true)) }


    //
    // MODULE: EXTRACT THE BUSCO GENES FOUND IN REFERENCE
    //
    GAWK_EXTRACT_BUSCOGENE (
        ch_busco_full_table,
        file("${projectDir}/bin/get_busco_gene.awk"),
        false
    )
    ch_versions                 = ch_versions.mix( GAWK_EXTRACT_BUSCOGENE.out.versions )


    //
    // LOGIC: ADDING LINE COUNT TO THE FILE FOR BETTER RESOURCE USAGE
    //
    GAWK_EXTRACT_BUSCOGENE.out.output
        .map { meta, file ->
            tuple ( [   id:     meta.id,
                        lines:  file.countLines()
                    ],
                    file
            )
        }
        .set { bedtools_input }
    //
    // MODULE: SORT THE EXTRACTED BUSCO GENE
    //
    BEDTOOLS_SORT(
        bedtools_input,
        []
    )
    ch_versions                 = ch_versions.mix( BEDTOOLS_SORT.out.versions )

    //
    // MODULE: CONVERT THE BED TO BIGBED
    //
    UCSC_BEDTOBIGBED(
        BEDTOOLS_SORT.out.sorted,
        dot_genome.map{it[1]},      // Gets file from tuple (meta, file)
        buscogene_as
    )
    ch_versions                 = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )

    //
    // SUBWORKFLOW: RUN ANCESTRAL BUSCO ID (ONLY AVAILABLE FOR LEPIDOPTERA)
    // LOGIC: AGGREGATE DATA AND FILTER ON CLASS
    //
    lineageinfo
        .combine(ch_busco_full_table)
        .combine(ancestral_table)
        .filter { lineage, _meta, _btable, _atable ->
            lineage.split('_')[0] == "lepidoptera"
        }
        .multiMap { _lineage, meta, busco_full_table, ancestral_table_ ->
            busco_table: tuple( meta, busco_full_table )
            atable:      ancestral_table_
        }
        .set{ ch_busco_lep_data }


    ANCESTRAL_GENE (
        ch_busco_lep_data.busco_table,
        dot_genome,
        buscogene_as,
        ch_busco_lep_data.atable
    )
    ch_versions                 = ch_versions.mix( ANCESTRAL_GENE.out.versions )

    emit:
    ch_buscogene_bigbed         = UCSC_BEDTOBIGBED.out.bigbed
    ch_ancestral_bigbed         = ANCESTRAL_GENE.out.ch_ancestral_bigbed
    versions                    = ch_versions

}
