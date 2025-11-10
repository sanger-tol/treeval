#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { FIND_TELOMERE_REGIONS         } from '../../../modules/local/find/telomere_regions/main'
include { GAWK as GAWK_SPLIT_DIRECTIONS } from '../../../modules/nf-core/gawk/main'

include { TELO_EXTRACTION               } from '../../../subworkflows/local/telo_extraction/main'

workflow TELO_FINDER {

    take:
    reference_tuple     // Channel: tuple [ val(meta), path(fasta) ]
    teloseq

    main:
    ch_versions     = Channel.empty()


    //
    // MODULE: FINDS THE TELOMERIC SEQEUNCE IN REFERENCE
    //
    FIND_TELOMERE_REGIONS (
        reference_tuple,
        teloseq
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_REGIONS.out.versions )

    FIND_TELOMERE_REGIONS.out.telomere
        .map{ meta, file ->
            def new_meta = meta + [direction: 0]
            [new_meta, file]
        }
        .set { ch_full_telomere }

    //
    // MODULE: SPLIT THE TELOMERE FILE INTO 5' and 3' FILES
    //              THIS IS RUNNING ON A LOCAL VERSION OF THE GAWK MODULE
    //
    if (params.split_telomere) {
        GAWK_SPLIT_DIRECTIONS (
            ch_full_telomere,
            [],
            true
        )
        ch_versions     = ch_versions.mix( GAWK_SPLIT_DIRECTIONS.out.versions )

        //
        // LOGIC: COLLECT FILES AND ITERATE THROUGH
        //          ADD DIRECTION BASED ON:
        //              0: FULL TELOMERE FILE
        //              3: FOR 3Prime DIRECTION
        //              5: For 5Prime DIRECTION
        //          THIS PRODUCES A TRIO OF CHANNELS: [meta], file
        //          FILTER FOR SIZE > 0 FOR SAFETY
        //
        GAWK_SPLIT_DIRECTIONS.out.output
            .flatMap { meta, files ->
                files
                    .findAll { file -> file.size() > 0 }
                    .collect { file ->
                        if (file.name.contains("direction.0")) {
                            new_meta = meta + [direction: 5]
                        }
                        if (file.name.contains("direction.1")) {
                            new_meta = meta + [direction: 3]
                        }
                        [new_meta, file]
                    }
            }
            .mix( ch_full_telomere )
            .set { ch_regions_for_extraction }


    } else {
        ch_regions_for_extraction  = ch_full_telomere
    }


    //
    // SUBWORKFLOW: TELO_EXTRACTION
    //              - The prime5.mix(prime3) creates a queue channel to execute
    //                  TELO_EXTRACTION per item in channel
    //
    TELO_EXTRACTION (
        ch_regions_for_extraction
    )
    ch_versions     = ch_versions.mix( TELO_EXTRACTION.out.versions )


    TELO_EXTRACTION.out.bedgraph_file
        .map{ _meta, bedgraph ->
            bedgraph
        }
        .collect()
        .set { telo_bedgraphs }


    emit:
    bed_file        = TELO_EXTRACTION.out.bed_file.collect()    // Not used anymore
    bed_gz_tbi      = TELO_EXTRACTION.out.bed_gz_tbi.collect()  // Not used anymore
    bedgraph_file   = telo_bedgraphs                            // Used in pretext_graph
    versions        = ch_versions
}
