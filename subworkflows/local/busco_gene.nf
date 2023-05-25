#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { BUSCO                          } from '../../modules/nf-core/busco/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UCSC_BEDTOBIGBED               } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                  } from '../../modules/nf-core/bedtools/sort/main'
include { EXTRACT_BUSCOGENE              } from '../../modules/local/extract_buscogene'

workflow BUSCO_GENE {
    take:
        reference_tuple      // channel: [id: sample_id], reference_file
        lineageinfo          // channel: [ meta, lineage_db, lineage ] 
        lineagespath         // channel: [ meta, /path/to/buscoDB, lineage ] 
        dot_genome           // Channel: [val(meta), [ datafile ]]
        buscogene_as         // channel val(dot_as location)

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

    BEDTOOLS_SORT(EXTRACT_BUSCOGENE.out.genefile, [])
    ch_versions             = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    UCSC_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, dot_genome.map{it[1]}, buscogene_as)
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_buscogene_bigbed    = UCSC_BEDTOBIGBED.out.bigbed
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
