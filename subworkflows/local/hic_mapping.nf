#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, pretext, .hic

nextflow.enable.dsl=2

// MODULE IMPORT
include { GENERATE_CRAM_CSV                                         } from '../../modules/local/generate_cram_csv'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT                    } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { BWAMEM2_INDEX                                             } from '../../modules/nf-core/bwamem2/index/main'
workflow HIC_MAPPING {
    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    hic_reads

    main:
    ch_versions         = Channel.empty()

    // bwamem2 indexing on reference output the folder of indexing files
    BWAMEM2_INDEX (reference_tuple)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    // make channel of hic reads as input of GENERATE_CRAM_CSV
    reference_tuple
        .combine( hic_reads )
        .map { meta, ref, hic_reads ->
                tuple([ id: meta.id, single_end: true], hic_reads) }
        .set { get_reads_input }

    ch_grab  = GrabFiles(get_reads_input)

    // generate a cram csv file contains required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    GENERATE_CRAM_CSV ( ch_grab )
    ch_versions = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)

    // organise all parametres into a channel for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    ch_filtering_input          = GENERATE_CRAM_CSV.out.csv
                            .splitCsv()
                            .combine (reference_tuple)
                            .combine (BWAMEM2_INDEX.out.index)
                            .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
                                  tuple([ id: cram_id.id], file(cram_info[0]), cram_info[1], file(ref_dir), cram_info[2], cram_info[3], cram_info[4], cram_info[5], cram_info[6], bwa_path.toString()+'/'+ref_dir.toString().split('/')[-1])
                                      
                            }
                            
    // parallel proccessing bwa-mem2 alignment by given interval of containers from cram files
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT ( ch_filtering_input  )
    ch_versions = ch_versions.mix(CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions)


    emit:
    mappedbam           = CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam
    versions            = ch_versions.ifEmpty(null)
}


process GrabFiles {

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*cram")

    "true"
}