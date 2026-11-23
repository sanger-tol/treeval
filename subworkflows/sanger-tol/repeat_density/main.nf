#!/usr/bin/env nextflow

//
// SANGER_TOL SUBWORKFLOWMODULE IMPORT BLOCK
//
include { REPEAT_MASKING                    } from '../repeat_masking/main'
include { FEATURE_DENSITY                   } from '../feature_density/main'

//
// MODULE IMPORT BLOCK
//
include { GAWK as GAWK_EXTRACT_REPEATS      } from '../../../modules/nf-core/gawk/main'


workflow REPEAT_DENSITY {
    take:
    ch_reference     // Channel: tuple [ val(meta), path(file) ]
    ch_chrom_sizes

    main:
    //
    // MODULE: MARK UP THE REPEAT REGIONS OF THE REFERENCE GENOME
    //
    REPEAT_MASKING (
        ch_reference
    )


    //
    // MODULE: USE USTAT OUTPUT TO EXTRACT REPEATS FROM FASTA
    //
    ch_extract_repeats_awk = channel.of('''\
        BEGIN { FS = " - "; OFS = "\\t" }
        /^>/ {
            header = substr($0, 2)
            next
        }
        {
            print header, $1, $2
        }'''.stripIndent())
        .collectFile(name: "extract_repeats.awk", cache: true)
        .collect()

    GAWK_EXTRACT_REPEATS(
        REPEAT_MASKING.out.repeat_intervals.map { meta, files -> tuple(meta, files, 'bed') },
        ch_extract_repeats_awk,
        false
    )


    //
    // SUBWORKFLOW: FINDE FEATURE ACROSS GENOME AND OUTPUT A DENSITY FILE
    //              AND TABIX INDEX
    //
    FEATURE_DENSITY(
        GAWK_EXTRACT_REPEATS.out.output,
        ch_chrom_sizes
    )

    emit:
    repeat_density       = FEATURE_DENSITY.out.density_file
    repeat_density_tabix = FEATURE_DENSITY.out.density_tabix
}
