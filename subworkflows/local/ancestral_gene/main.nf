#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { EXTRACT_ANCESTRAL  } from '../../../modules/local/extract/ancestral/main'
include { ASSIGN_ANCESTRAL   } from '../../../modules/local/assign/ancestral/main'
include { BEDTOOLS_SORT      } from '../../../modules/nf-core/bedtools/sort/main'
include { UCSC_BEDTOBIGBED   } from '../../../modules/nf-core/ucsc/bedtobigbed/main'

workflow ANCESTRAL_GENE {
    take:
    busco_full_table     // Channel: tuple [val(meta),/path/to/busco/output/**/fulltable.tsv ]
    dot_genome           // Channel: tuple [val(meta), [ datafile ]]
    buscogene_as         // Channel: val(dot_as location)
    ancestral_table      // Channel: val(ancestral_table location)

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: EXTRACTS ANCESTRALLY LINKED BUSCO GENES FROM FULL TABLE
    //
    EXTRACT_ANCESTRAL(
        busco_full_table,
        ancestral_table
    )
    ch_versions             = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)

    //
    // MODULE: ASSIGN EXTRACTED GENES TO ANCESTRAL GROUPS
    //
    ASSIGN_ANCESTRAL(
        EXTRACT_ANCESTRAL.out.comp_location,
        busco_full_table.map { _meta, fulltable -> fulltable }
    )
    ch_versions             = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)

    //
    // MODULES: SORT THE BED FILE
    //
    BEDTOOLS_SORT(
        ASSIGN_ANCESTRAL.out.assigned_bed,
        []
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    //
    // MODULES: CONVERT BED TO INDEXED BIGBED
    //
    UCSC_BEDTOBIGBED(
        BEDTOOLS_SORT.out.sorted,
        dot_genome.map{ it[1] },      // Pull file from tuple(meta, file)
        buscogene_as
    )
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_ancestral_bigbed     = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions
}
