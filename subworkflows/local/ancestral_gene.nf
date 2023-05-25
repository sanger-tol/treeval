#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { EXTRACT_ANCESTRAL  } from '../../modules/local/extract_ancestral'
include { ASSIGN_ANCESTRAL   } from '../../modules/local/assign_ancestral'
include { BEDTOOLS_SORT      } from '../../modules/nf-core/bedtools/sort/main'
include { UCSC_BEDTOBIGBED   } from '../../modules/nf-core/ucsc/bedtobigbed/main'

workflow ANCESTRAL_GENE {
    take:
        busco_dir            // Channel: [val(meta),/path/to/busco/output/dir]
        dot_genome           // Channel: [val(meta), [ datafile ]]
        buscogene_as         // channel val(dot_as location)
        ancestral_table      // channel val(ancestral_table location)

    main:
    ch_versions             = Channel.empty()

    ch_grab  = GrabFiles(busco_dir)

    EXTRACT_ANCESTRAL(ch_grab, ancestral_table)
    ch_versions = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)

    ch_grab
        .map { meta, fulltable
                -> fulltable
            }
        .set { assignanc_input }

    ASSIGN_ANCESTRAL(EXTRACT_ANCESTRAL.out.comp_location, assignanc_input )
    ch_versions = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)

    BEDTOOLS_SORT(ASSIGN_ANCESTRAL.out.assigned_bed, [])
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    UCSC_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, dot_genome.map{it[1]}, buscogene_as)
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_ancestral_bigbed    = UCSC_BEDTOBIGBED.out.bigbed
    versions               = ch_versions.ifEmpty(null)
}
process GrabFiles {

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*/*/full_table.tsv")

    "true"
}