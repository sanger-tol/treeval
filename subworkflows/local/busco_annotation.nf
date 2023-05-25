#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { BUSCO_GENE } from './busco_gene'
include { ANCESTRAL_GENE } from './ancestral_gene'

workflow BUSCO_ANNOTATION {
    take:
    dot_genome           // channel: [val(meta), [ datafile ]]
    reference_tuple      // channel: [val(meta), [ datafile ]]
    assembly_classT      // channel: val(class)
    lineageinfo          // channel: val(lineage_db)
    lineagespath         // channel: val(/path/to/buscoDB)
    buscogene_as         // channel: val(dot_as location)
    ancestral_table      // channel: val(ancestral_table location)

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: TAKES A SINGLE LIKE CSV STRING AND CONVERTS TO LIST OF VALUES
    //          LIST IS MERGED WITH DATA_DIRECTORY AND ORGANISM_CLASS
    //

    reference_tuple
            .combine( lineageinfo )
            .combine( lineagespath )
            .combine( assembly_classT )
            .combine( ancestral_table )
            .branch {
                lep:     it[4] == "lepidoptera"
                general: it[4] != "lepidoptera"
            }
            .set{ch_busco_data}

    ch_busco_data.lep
            .multiMap { data ->
                    refT:   tuple(data[0], data[1])
                    liT:    data[2]
                    lpT:    data[3]
                    atable: data[5]
            }
            .set{ch_busco_lep_data}

    ch_busco_data.general
            .multiMap { data ->
                    refT:   tuple(data[0], data[1])
                    liT:    data[2]
                    lpT:    data[3]
            }
            .set{ch_busco_general_data}
    
    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR RNA, NUCLEAR AND COMPLEMENT_DNA DATA, EMITS BIGBED
    //
    BUSCO_GENE (ch_busco_general_data.refT,
                ch_busco_general_data.liT,
                ch_busco_general_data.lpT,
                dot_genome,
                buscogene_as )
    ch_versions = ch_versions.mix(BUSCO_GENE.out.versions)

    ANCESTRAL_GENE (ch_busco_lep_data.refT,
                    ch_busco_lep_data.liT,
                    ch_busco_lep_data.lpT,
                    dot_genome,
                    buscogene_as,
                    ch_busco_lep_data.atable)
    ch_versions = ch_versions.mix(ANCESTRAL_GENE.out.versions)

    emit:
    ch_buscogene_bigbed        = BUSCO_GENE.out.bigbed
    ch_ancestral_busco_bigbed  = ANCESTRAL_GENE.out.ch_busco_bigbed
    ch_ancestral_bigbed        = ANCESTRAL_GENE.out.ch_ancestral_bigbed
    versions = ch_versions

}
