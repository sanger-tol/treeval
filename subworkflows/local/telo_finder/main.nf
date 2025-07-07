#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { GAWK as GAWK_UPPER_SEQUENCE   } from '../../../modules/nf-core/gawk/main'
include { FIND_TELOMERE_REGIONS         } from '../../../modules/local/find/telomere_regions/main'
include { GAWK as GAWK_SPLIT_DIRECTIONS } from '../../../modules/local/gawk/main'

include { TELO_EXTRACTION               } from '../../../subworkflows/local/telo_extraction/main'

workflow TELO_FINDER {

    take:
    reference_tuple     // Channel: tuple [ val(meta), path(fasta) ]
    teloseq

    main:
    ch_versions     = Channel.empty()

    //
    // MODULE: UPPERCASE THE REFERENCE SEQUENCE
    //
    GAWK_UPPER_SEQUENCE(
        reference_tuple,
        [],
        false,
    )
    ch_versions     = ch_versions.mix( GAWK_UPPER_SEQUENCE.out.versions )

    //
    // MODULE: FINDS THE TELOMERIC SEQEUNCE IN REFERENCE
    //
    FIND_TELOMERE_REGIONS (
        GAWK_UPPER_SEQUENCE.out.output,
        teloseq
    )
    ch_versions     = ch_versions.mix( FIND_TELOMERE_REGIONS.out.versions )


    //
    // MODULE: SPLIT THE TELOMERE FILE INTO 5' and 3' FILES
    //              THIS IS RUNNING ON A LOCAL VERSION OF THE GAWK MODULE
    //
    if (params.split_telomere) {
        GAWK_SPLIT_DIRECTIONS (
            FIND_TELOMERE_REGIONS.out.telomere,
            file("${projectDir}/bin/gawk_split_directions.awk"),
            false
        )
        ch_versions     = ch_versions.mix( GAWK_SPLIT_DIRECTIONS.out.versions )

        GAWK_SPLIT_DIRECTIONS.out.prime5
            .map { meta, file ->
                 println file
                 def prime = file.name.contains(".0.") ? "5P" : "NA"
                 tuple( [id: meta.id + "_" + prime], file)
            }
            .set { prime5_telo }

        GAWK_SPLIT_DIRECTIONS.out.prime3
            .map { meta, file ->
                println file
                def prime = file.name.contains(".1.") ? "3P" : "NA"
                tuple( [id: meta.id + "_" + prime], file)
            }
            .set { prime3_telo }

        telo_for_extraction = prime5_telo.mix(prime3_telo)
        telo_for_extraction.view{"SPLIT TELOS: $it"}

    } else {
        telo_for_extraction = FIND_TELOMERE_REGIONS.out.telomere
    }

    //
    // SUBWORKFLOW: TELO_EXTRACTION
    //              - The prime5.mix(prime3) creates a queue channel to execute
    //                  TELO_EXTRACTION per item in channel
    //
    TELO_EXTRACTION (
        telo_for_extraction
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
