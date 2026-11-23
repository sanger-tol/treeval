#!/usr/bin/env nextflow

//
// SANGER_TOL SUBWORKFLOW IMPORT BLOCK
//
include { GAP_FINDER                        } from '../gap_finder/main'
include { TELO_FINDER                       } from '../telo_finder/main'
include { READ_COVERAGE                     } from '../read_coverage/main'
include { REPEAT_DENSITY                    } from '../repeat_density/main'

workflow PRETEXT_ACCESSORY_FILES {
    take:
    ch_reference_tuple  // Channel [ val(meta), path(file)   ]
    ch_reference_sizes  // Channel [ val(meta), path(file)   ]
    ch_longread_reads   // Channel [ val(meta), [ path(file) ] ]
    ch_teloseq          // Channel [ val(meta), val(telomere_sequence) ]
    val_run_telomere    // val(bool)
    val_run_gaps        // val(bool)
    val_run_coverage    // val(bool)
    val_run_repeat_den  // val(bool)
    _val_run_busco_track // val(bool) PLACEHOLDER
    _val_run_pebble      // val(bool) PLACEHOLDER
    val_output_indexes  // val(bool)


    main:

    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    GAP_FINDER(
        ch_reference_tuple.filter { val_run_gaps },
        val_output_indexes
    )


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH LONGREAD READS AND REFERENCE
    //
    TELO_FINDER(
        ch_reference_tuple.filter { val_run_telomere },
        ch_teloseq,
        val_output_indexes
    )

    telomere_data = TELO_FINDER.out.telomere_bed
        .mix( TELO_FINDER.out.telomere_bed_fwd )
        .mix( TELO_FINDER.out.telomere_bed_rev )

    telomere_windows = TELO_FINDER.out.windows_all
        .mix( TELO_FINDER.out.windows_fwd )
        .mix( TELO_FINDER.out.windows_rev )


    //
    // LOGIC: MAP TOGETHER THE REF, SIZES AND READS SO WE DON'T GET MISMATCHED INPUTS
    //
    data = ch_reference_tuple
        .combine(ch_reference_sizes
            .map{ meta, file -> tuple([id: meta.id], file)}, by: 0)
        .combine(ch_longread_reads
            .map{ meta, file_list -> tuple([id: meta.id], file_list)}, by: 0)
        .multiMap { meta, ref, sizes, reads_list ->
            reference: tuple(meta, ref)
            sizes: tuple(meta, sizes)
            longreads: tuple(meta, reads_list)
        }


    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    REPEAT_DENSITY(
        data.reference.filter { val_run_repeat_den },
        data.sizes
    )


    //
    // SUBWORKFLOW: Takes reference, longread reads
    //
    READ_COVERAGE(
        data.longreads,
        data.reference.filter { val_run_coverage },
        data.sizes.map { _meta, file -> file }
    )


    //
    // MODULE: PEBBLE PLACEHOLDER
    //
    // PEBBLE (
    //  ch_reference_tuple
    // )


    //
    // SUBWORKFLOW: BUSCO_TRACK PLACEHOLDER
    //
    // BUSCO_TRACK (
    //  ch_reference_tuple,
    //  busco_data
    //)


    emit:
    gap_file            = GAP_FINDER.out.gap_file
    gap_bed_and_index   = GAP_FINDER.out.gap_and_index
    repeat_file         = REPEAT_DENSITY.out.repeat_density
    telo_file           = telomere_data                     // This is the possible collection of telomere files
    telo_windows        = telomere_windows                  // This is the possible collection of window files
    telo_gz_and_index   = TELO_FINDER.out.gz_index
    coverage_file       = READ_COVERAGE.out.bigwig
}
