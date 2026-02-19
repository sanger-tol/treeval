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
include { UCSC_BEDGRAPHTOBIGWIG             } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GNU_SORT as GNU_SORT_A            } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_B            } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_C            } from '../../../modules/nf-core/gnu/sort/main'
include { GAWK as GAWK_RENAME_IDS           } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REPLACE_DOTS         } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REFORMAT_INTERSECT   } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIPTABIX                  } from '../../../modules/nf-core/tabix/bgziptabix'


workflow REPEAT_DENSITY {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    dot_genome

    main:
    ch_versions         = channel.empty()
    //
    // MODULE: MARK UP THE REPEAT REGIONS OF THE REFERENCE GENOME
    //
    WINDOWMASKER_MKCOUNTS (
        reference_tuple
    )

    //
    // MODULE: CALCULATE THE STATISTICS OF THE MARKED UP REGIONS
    //
    WINDOWMASKER_USTAT(
        WINDOWMASKER_MKCOUNTS.out.counts,
        reference_tuple
    )

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

    //
    // LOGIC: COMBINE TWO CHANNELS AND OUTPUT tuple(meta, windows_file, repeat_file)
    //
    BEDTOOLS_MAKEWINDOWS.out.bed
        .combine( EXTRACT_REPEAT.out.bed )
        .map{ meta, windows_file, _repeat_meta, repeat_file ->
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

    //
    // MODULE: FIXES IDS FOR REPEATS
    //
    GAWK_RENAME_IDS(
        BEDTOOLS_INTERSECT.out.intersect,
        file("${projectDir}/bin/gawk_replace_dots.awk"),
        false
    )

    //
    // MODULE: SORTS THE ABOVE BED FILES
    //
    GNU_SORT_A (
        GAWK_RENAME_IDS.out.output      // Intersect file
    )

    dot_genome                          // Rename file so input/output are different
        .map { meta, sizes ->
            tuple(
                meta,
                file(sizes).moveTo("${meta.id}.unsorted.genome")
                )
        }
        .set { ch_sizes }        

    GNU_SORT_B (
        ch_sizes                      // Genome file - Will not run unless genome file is sorted to
    )

    GNU_SORT_C (
        BEDTOOLS_MAKEWINDOWS.out.bed    // Windows file
    )

    //
    // MODULE: ADDS 4TH COLUMN TO BED FILE USED IN THE REPEAT DENSITY GRAPH
    //
    GAWK_REFORMAT_INTERSECT (
        GNU_SORT_A.out.sorted,
        file("${projectDir}/bin/gawk_reformat_intersect.awk"),
        false
    )

    //
    // MODULE: TABIX AND GZIP THE REPEAT DENSITY BED FILE FOR JBROWSE
    //
    TABIX_BGZIPTABIX (
        GAWK_REFORMAT_INTERSECT.out.output
    )

    //
    // LOGIC: COMBINES THE REFORMATTED INTERSECT FILE AND WINDOWS FILE CHANNELS AND SORTS INTO
    //        tuple(intersect_meta, windows file, intersect file)
    //
    GAWK_REFORMAT_INTERSECT.out.output
        .combine( GNU_SORT_C.out.sorted )
        .map{ intersect_meta, bed, _sorted_meta, windows_file ->
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

    //
    // MODULE: REPLACES . WITH 0 IN MAPPED FILE
    //
    GAWK_REPLACE_DOTS (
        BEDTOOLS_MAP.out.mapped,
        file("${projectDir}/bin/gawk_replace_dots.awk"),
        false
    )

    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGWIG FILE
    //
    UCSC_BEDGRAPHTOBIGWIG(
        GAWK_REPLACE_DOTS.out.output,
        GNU_SORT_B.out.sorted.map { _meta, file -> file } // Pulls file from tuple of meta and file
    )

    emit:
    repeat_density      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bed_gz_tbi          = TABIX_BGZIPTABIX.out.gz_index
    versions            = ch_versions
}
