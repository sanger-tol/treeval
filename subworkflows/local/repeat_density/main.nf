#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { WINDOWMASKER_USTAT                } from '../../../modules/nf-core/windowmasker/ustat/main'
include { WINDOWMASKER_MKCOUNTS             } from '../../../modules/nf-core/windowmasker/mkcounts/main'
include { EXTRACT_REPEAT                    } from '../../../modules/local/extract/repeat/main'
include { BEDTOOLS_INTERSECT                } from '../../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_MAKEWINDOWS              } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAP                      } from '../../../modules/nf-core/bedtools/map/main'
include { RENAME_IDS                        } from '../../../modules/local/rename/ids/main'
include { UCSC_BEDGRAPHTOBIGWIG             } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GNU_SORT as GNU_SORT_A            } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_B            } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_C            } from '../../../modules/nf-core/gnu/sort/main'
include { REFORMAT_INTERSECT                } from '../../../modules/local/reformat/intersect/main'
include { REPLACE_DOTS                      } from '../../../modules/local/replace/dots/main'
include { TABIX_BGZIPTABIX                  } from '../../../modules/nf-core/tabix/bgziptabix'

workflow REPEAT_DENSITY {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    dot_genome

    main:
    ch_versions         = Channel.empty()
    //
    // MODULE: MARK UP THE REPEAT REGIONS OF THE REFERENCE GENOME
    //
    WINDOWMASKER_MKCOUNTS (
        reference_tuple
    )
    ch_versions         = ch_versions.mix( WINDOWMASKER_MKCOUNTS.out.versions )

    //
    // MODULE: CALCULATE THE STATISTICS OF THE MARKED UP REGIONS
    //
    WINDOWMASKER_USTAT(
        WINDOWMASKER_MKCOUNTS.out.counts,
        reference_tuple
    )
    ch_versions         = ch_versions.mix( WINDOWMASKER_USTAT.out.versions )

    //
    // MODULE: USE USTAT OUTPUT TO EXTRACT REPEATS FROM FASTA
    //
    EXTRACT_REPEAT(
        WINDOWMASKER_USTAT.out.intervals
    )
    ch_versions         = ch_versions.mix( EXTRACT_REPEAT.out.versions )

    //
    // MODULE: CREATE WINDOWS FROM .GENOME FILE
    //
    BEDTOOLS_MAKEWINDOWS(
        dot_genome
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_MAKEWINDOWS.out.versions )

    //
    // LOGIC: COMBINE TWO CHANNELS AND OUTPUT tuple(meta, windows_file, repeat_file)
    //
    BEDTOOLS_MAKEWINDOWS.out.bed
        .combine( EXTRACT_REPEAT.out.bed )
        .map{ meta, windows_file, repeat_meta, repeat_file ->
                    tuple (
                        meta,
                        windows_file,
                        repeat_file
                    )
        }
        .set { intervals }

    //
    // MODULE: GENERATES THE REPEAT FILE FROM THE WINDOW FILE AND GENOME FILE
    //
    BEDTOOLS_INTERSECT(
        intervals,
        dot_genome
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_INTERSECT.out.versions )

    //
    // MODULE: FIXES IDS FOR REPEATS
    //
    RENAME_IDS(
        BEDTOOLS_INTERSECT.out.intersect
    )
    ch_versions         = ch_versions.mix( RENAME_IDS.out.versions )

    //
    // MODULE: SORTS THE ABOVE BED FILES
    //
    GNU_SORT_A (
        RENAME_IDS.out.bed              // Intersect file
    )
    ch_versions         = ch_versions.mix( GNU_SORT_A.out.versions )

    GNU_SORT_B (
        dot_genome                      // Genome file - Will not run unless genome file is sorted to
    )
    ch_versions         = ch_versions.mix( GNU_SORT_B.out.versions )

    GNU_SORT_C (
        BEDTOOLS_MAKEWINDOWS.out.bed    // Windows file
    )
    ch_versions         = ch_versions.mix( GNU_SORT_C.out.versions )

    //
    // MODULE: ADDS 4TH COLUMN TO BED FILE USED IN THE REPEAT DENSITY GRAPH
    //
    REFORMAT_INTERSECT (
        GNU_SORT_A.out.sorted
    )
    ch_versions         = ch_versions.mix( GNU_SORT_C.out.versions )

    //
    // MODULE: TABIX AND GZIP THE REPEAT DENSITY BED FILE FOR JBROWSE
    //
    TABIX_BGZIPTABIX (
        REFORMAT_INTERSECT.out.bed
    )
    ch_versions     = ch_versions.mix( TABIX_BGZIPTABIX.out.versions )

    //
    // LOGIC: COMBINES THE REFORMATTED INTERSECT FILE AND WINDOWS FILE CHANNELS AND SORTS INTO
    //        tuple(intersect_meta, windows file, intersect file)
    //
    REFORMAT_INTERSECT.out.bed
        .combine( GNU_SORT_C.out.sorted )
        .map{ intersect_meta, bed, sorted_meta, windows_file ->
                    tuple (
                        intersect_meta,
                        windows_file,
                        bed
                    )
        }
        .set { for_mapping }

    //
    // MODULE: MAPS THE REPEATS AGAINST THE REFERENCE GENOME
    //
    BEDTOOLS_MAP(
        for_mapping,
        GNU_SORT_B.out.sorted
    )
    ch_versions         = ch_versions.mix( BEDTOOLS_MAP.out.versions )

    //
    // MODULE: REPLACES . WITH 0 IN MAPPED FILE
    //
    REPLACE_DOTS (
        BEDTOOLS_MAP.out.mapped
    )
    ch_versions         = ch_versions.mix( REPLACE_DOTS.out.versions )

    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGWIG FILE
    //
    UCSC_BEDGRAPHTOBIGWIG(
        REPLACE_DOTS.out.bed,
        GNU_SORT_B.out.sorted.map { it[1] } // Pulls file from tuple of meta and file
    )
    ch_versions         = ch_versions.mix( UCSC_BEDGRAPHTOBIGWIG.out.versions )

    emit:
    repeat_density      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bed_gz_tbi          = TABIX_BGZIPTABIX.out.gz_tbi
    versions            = ch_versions.ifEmpty(null)
}
