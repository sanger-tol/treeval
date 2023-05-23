#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { BUSCO                                           } from '../../modules/nf-core/busco/main'
include { SAMTOOLS_FAIDX                                  } from '../../modules/nf-core/samtools/faidx/main'
include { UCSC_BEDTOBIGBED as UCSC_BEDTOBIGBED_BUSCO      } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_BUSCO            } from '../../modules/nf-core/bedtools/sort/main'
include { EXTRACT_BUSCOGENE                               } from '../../modules/local/extract_buscogene'
include { EXTRACT_ANCESTRAL                               } from '../../modules/local/extract_ancestral'
include { ASSIGN_ANCESTRAL                                } from '../../modules/local/assign_ancestral'
include { BEDTOOLS_SORT as BEDTOOLS_SORT_ANCESTRAL        } from '../../modules/nf-core/bedtools/sort/main'
include { UCSC_BEDTOBIGBED as UCSC_BEDTOBIGBED_ANCESTRAL  } from '../../modules/nf-core/ucsc/bedtobigbed/main'

workflow ANCESTRAL_GENE {
    take:
        reference_tuple      // channel: [id: sample_id], reference_file
        lineageinfo          // channel: [ meta, lineage_db, lineage ] 
        lineagespath         // channel: [ meta, /path/to/buscoDB, lineage ] 
        dot_genome           // Channel: [val(meta), [ datafile ]]
        buscogene_as         // channel val(dot_as location)
        ancestral_table      // channel val(ancestral_table location)
    


    main:
    ch_versions             = Channel.empty()
     
    // 
    // MODULE: RUN BUSCO TO OBTAIN FULL_TABLE.CSV
    //         EMITS FULL_TABLE.CSV
    //
    BUSCO ( reference_tuple, 
            lineageinfo,
            lineagespath,
            [] )
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())

    ch_grab  = GrabFiles(BUSCO.out.busco_dir)

    EXTRACT_BUSCOGENE (ch_grab)
    ch_versions = ch_versions.mix(EXTRACT_BUSCOGENE.out.versions)

    BEDTOOLS_SORT_BUSCO(EXTRACT_BUSCOGENE.out.genefile, [])
    ch_versions             = ch_versions.mix(BEDTOOLS_SORT_BUSCO.out.versions)

    UCSC_BEDTOBIGBED_BUSCO(BEDTOOLS_SORT_BUSCO.out.sorted, dot_genome.map{it[1]}, buscogene_as)
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED_BUSCO.out.versions)


    EXTRACT_ANCESTRAL(ch_grab, ancestral_table)
    ch_versions = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)


    ch_grab
        .map { meta, fulltable
                -> fulltable
            }
        .set { assignanc_input }

    ASSIGN_ANCESTRAL(EXTRACT_ANCESTRAL.out.comp_location, assignanc_input )
    ch_versions = ch_versions.mix(EXTRACT_ANCESTRAL.out.versions)

    BEDTOOLS_SORT_ANCESTRAL(ASSIGN_ANCESTRAL.out.assigned_bed, [])
    ch_versions = ch_versions.mix(BEDTOOLS_SORT_ANCESTRAL.out.versions)

    UCSC_BEDTOBIGBED_ANCESTRAL(BEDTOOLS_SORT_ANCESTRAL.out.sorted, dot_genome.map{it[1]}, buscogene_as)
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED_ANCESTRAL.out.versions)


    emit:
    ch_busco_bigbed        = UCSC_BEDTOBIGBED_BUSCO.out.bigbed
    ch_ancestral_bigbed    = UCSC_BEDTOBIGBED_ANCESTRAL.out.bigbed
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