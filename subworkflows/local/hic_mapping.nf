#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

//
// MODULE IMPORT BLOCK
//
include { COOLER_CLOAD                                    } from '../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY                                  } from '../../modules/nf-core/cooler/zoomify/main'
include { PRETEXTMAP as PRETEXTMAP_STANDRD                } from '../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_HIGHRES                } from '../../modules/nf-core/pretextmap/main'
include { PRETEXTSNAPSHOT as SNAPSHOT_SRES                } from '../../modules/nf-core/pretextsnapshot/main'
include { GENERATE_CRAM_CSV                               } from '../../modules/local/generate_cram_csv'
include { JUICER_TOOLS_PRE                                } from '../../modules/local/juicer_tools_pre'
include { SUBSAMPLE_BAM                                   } from '../../modules/local/subsample_bam.nf'
include { PRETEXT_INGESTION as PRETEXT_INGEST_SNDRD       } from '../../subworkflows/local/pretext_ingestion'
include { PRETEXT_INGESTION as PRETEXT_INGEST_HIRES       } from '../../subworkflows/local/pretext_ingestion'
include { HIC_BAMTOBED as HIC_BAMTOBED_COOLER             } from '../../subworkflows/local/hic_bamtobed'
include { HIC_BAMTOBED as HIC_BAMTOBED_JUICER             } from '../../subworkflows/local/hic_bamtobed'
include { HIC_MINIMAP2                                    } from '../../subworkflows/local/hic_minimap2'
include { HIC_BWAMEM2                                     } from '../../subworkflows/local/hic_bwamem2'

workflow HIC_MAPPING {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path( file )      ]
    reference_index     // Channel: tuple [ val(meta), path( file )      ]
    dot_genome          // Channel: tuple [ val(meta), path( datafile )  ]
    hic_reads_path      // Channel: tuple [ val(meta), path( directory ) ]
    assembly_id         // Channel: val( id )
    gap_file            // Channel: tuple [ val(meta), path( file )      ]
    coverage_file       // Channel: tuple [ val(meta), path( file )      ]
    avgcoverage_file    // Channel: tuple [ val(meta), path( file )      ]
    telo_file           // Channel: tuple [ val(meta), path( file )      ]
    repeat_density_file // Channel: tuple [ val(meta), path( file )      ]
    workflow_setting    // Channel: val( { RAPID | FULL | RAPID_TOL } )

    main:
    ch_versions         = Channel.empty()

    // COMMENT: 1000bp BIN SIZE INTERVALS FOR CLOAD
    ch_cool_bin         = Channel.of(1000)


    //
    // LOGIC: make channel of hic reads as input for GENERATE_CRAM_CSV
    //
    reference_tuple
        .combine(hic_reads_path)
        .map {meta, ref, hic_meta, hic_reads_path ->
                tuple(
                    [ id: meta.id, single_end: true],
                    hic_reads_path
                )
        }
        .set {get_reads_input}

    //
    // MODULE: generate a cram csv file containing the required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV (
        get_reads_input
    )
    ch_versions         = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)

    //
    // LOGIC: make branches for different hic aligner.
    //
    hic_reads_path
        .combine(reference_tuple)
        .map{meta, hic_read_path, ref_meta, ref->
             tuple(
                [   id : ref_meta,
                    aligner : meta.aligner
                ],
                ref
             )
        }
        .branch{
            minimap2      : it[0].aligner == "minimap2"
            bwamem2       : it[0].aligner == "bwamem2"
        }
        .set{ch_aligner}

    //
    // SUBWORKFLOW: mapping hic reads using minimap2
    //
    HIC_MINIMAP2 (
        ch_aligner.minimap2,
        GENERATE_CRAM_CSV.out.csv,
        reference_index
    )
    ch_versions         = ch_versions.mix(HIC_MINIMAP2.out.versions)
    mergedbam           = HIC_MINIMAP2.out.mergedbam

    //
    // SUBWORKFLOW: mapping hic reads using bwamem2
    //
    HIC_BWAMEM2 (
        ch_aligner.bwamem2,
        GENERATE_CRAM_CSV.out.csv,
        reference_index
    )
    ch_versions         = ch_versions.mix(HIC_BWAMEM2.out.versions)
    mergedbam           = mergedbam.mix(HIC_BWAMEM2.out.mergedbam)

    //
    // LOGIC: PREPARING PRETEXT MAP INPUT
    //
    mergedbam
        .combine(reference_tuple)
        .combine (dot_genome)
        .multiMap { bam_meta, bam, ref_meta, ref_fa, genome_meta, genome_file ->
            input_bam:  tuple( [    id: bam_meta.id,
                                    sz: file(bam).size() ],
                                bam
                        )
            // NOTE: Inject the genome file into the channel to speed up PretextMap
            reference:  tuple(  ref_meta,
                                ref_fa,
                                genome_file
                        )
        }
        .set {pretext_input}

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR LOW RES
    //
    PRETEXTMAP_STANDRD (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix(PRETEXTMAP_STANDRD.out.versions)

    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //
    PRETEXT_INGEST_SNDRD (
        PRETEXTMAP_STANDRD.out.pretext,
        gap_file,
        coverage_file,
        avgcoverage_file,
        telo_file,
        repeat_density_file
    )
    ch_versions         = ch_versions.mix(PRETEXT_INGEST_SNDRD.out.versions)

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
    //
    PRETEXTMAP_HIGHRES (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix(PRETEXTMAP_HIGHRES.out.versions)

    //
    // NOTICE: This could fail on LARGE hires maps due to some memory parameter in the C code
    //         of pretext graph. There is a "fixed" version in sanger /software which may need
    //         to be released in this case
    //
    PRETEXT_INGEST_HIRES (
        PRETEXTMAP_HIGHRES.out.pretext,
        gap_file,
        coverage_file,
        avgcoverage_file,
        telo_file,
        repeat_density_file
    )
    ch_versions         = ch_versions.mix(PRETEXT_INGEST_HIRES.out.versions)

    //
    // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    //
    SNAPSHOT_SRES (
        PRETEXTMAP_STANDRD.out.pretext
    )
    ch_versions         = ch_versions.mix (SNAPSHOT_SRES.out.versions)

    //
    // LOGIC: BRANCH TO SUBSAMPLE BAM IF LARGER THAN 50G
    //
    mergedbam
        .map{ meta, bam ->
            tuple(
                [   id : meta.id,
                    sz : file(bam).size()
                ],
                bam
            )
        }
        .branch {
            tosubsample    : it[0].sz >= 50000000000
            unmodified     : it[0].sz < 50000000000
        }
        .set {ch_merged_bam}

    // LOGIC: PREPARE BAMTOBED JUICER INPUT.
    if (workflow_setting != "RAPID_TOL" && params.juicer == false) {
        //
        // LOGIC: BRANCH TO SUBSAMPLE BAM IF LARGER THAN 50G
        //
       mergedbam
            .map{ meta, bam ->
                tuple(
                        [   id : meta.id,
                        sz : file(bam).size()
                    ],
                    bam
                )
            }
            .branch {
                tosubsample    : it[0].sz >= 50000000000
                unmodified     : it[0].sz < 50000000000
            }
                .set {ch_merged_bam}

        //
        // MODULE: SUBSAMPLE BAM
        //
        SUBSAMPLE_BAM (
            ch_merged_bam.tosubsample
        )
        ch_versions = ch_versions.mix (SUBSAMPLE_BAM.out.versions)

        //
        // LOGIC: COMBINE BRANCHED TO SINGLE OUTPUT
        //
        ch_subsampled_bam = SUBSAMPLE_BAM.out.subsampled_bam
        ch_subsampled_bam.mix(ch_merged_bam.unmodified)

        //
        // LOGIC: PREPARE BAMTOBED JUICER INPUT
        //
        ch_subsampled_bam
            .combine(reference_tuple)
            .multiMap {  meta, subsampled_bam, meta_ref, ref ->
                bam            :   tuple(meta, subsampled_bam )
                reference      :   tuple(meta_ref, ref)
            }
            .set {ch_bamtobed_juicer_input}

        //
        // SUBWORKFLOW: BAM TO BED FOR JUICER - USES THE SUBSAMPLED MERGED BAM
        //
        HIC_BAMTOBED_JUICER(
            ch_bamtobed_juicer_input.bam,
            ch_bamtobed_juicer_input.reference
        )
        ch_versions         = ch_versions.mix(HIC_BAMTOBED_JUICER.out.versions)

        //
        // LOGIC: PREPARE JUICER TOOLS INPUT
        //
        HIC_BAMTOBED_JUICER.out.paired_contacts_bed
            .combine(dot_genome)
            .multiMap {  meta, paired_contacts, meta_my_genome, my_genome ->
                paired      :   tuple([id: meta.id, single_end: true], paired_contacts)
                genome      :   my_genome
                id          :   meta.id
            }
            .set {ch_juicer_input}

        //
        // MODULE: GENERATE HIC MAP, ONLY IS PIPELINE IS RUNNING ON ENTRY FULL
        //
        JUICER_TOOLS_PRE(
            ch_juicer_input.paired,
            ch_juicer_input.genome,
            ch_juicer_input.id
        )
        ch_versions         = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)
    }

    //
    // LOGIC: PREPARE BAMTOBED COOLER INPUT
    //
    mergedbam
        .combine(reference_tuple)
        .multiMap {  meta, merged_bam, meta_ref, ref ->
            bam            :   tuple(meta, merged_bam )
            reference      :   tuple(meta_ref, ref)
        }
        .set {ch_bamtobed_cooler_input}

    //
    // SUBWORKFLOW: BAM TO BED FOR COOLER
    //
    HIC_BAMTOBED_COOLER(
        ch_bamtobed_cooler_input.bam,
        ch_bamtobed_cooler_input.reference
    )
    ch_versions         = ch_versions.mix(HIC_BAMTOBED_COOLER.out.versions)

    //
    // LOGIC: BIN CONTACT PAIRS
    //
    HIC_BAMTOBED_COOLER.out.paired_contacts_bed
        .join(HIC_BAMTOBED_COOLER.out.sorted_bed)
        .combine( h_cool_bin)
        .set {ch_binned_pairs}

    //
    // LOGIC: PREPARE COOLER INPUT
    //
    ch_binned_pairs
        .combine(dot_genome)
        .multiMap {meta, pairs, bed, cool_bin, meta_my_genome, my_genome ->
            cooler_in   : tuple (meta, pairs, bed, cool_bin)
            genome_file : my_genome
        }
        .set {ch_cooler}

    //
    // MODULE: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    //
    COOLER_CLOAD(
        ch_cooler.cooler_in,
        ch_cooler.genome_file
    )
    ch_versions         = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // LOGIC: REFACTOR CHANNEL FOR ZOOMIFY
    //
    COOLER_CLOAD.out.cool
        .map{meta, cools, cool_bin ->
            [meta, cools]
        }
        .set{ch_cool}

    //
    // MODULE: ZOOM COOL TO MCOOL
    //
    COOLER_ZOOMIFY(ch_cool)
    ch_versions         = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    //
    // LOGIC: FOR REPORTING
    //

    ch_cram_files = GrabFiles( hic_reads_path )

    ch_cram_files
        .collect()
        .map {meta, cram ->
            tuple( [    id: 'cram',
                        sz: cram instanceof ArrayList ? cram.collect { it.size()} : cram.size(),
                    ],
                    cram
            )
        }
        .combine(GENERATE_CRAM_CSV.out.csv)
        .map { meta, data, meta2, csv ->
            tuple( [    id: meta.id,
                        sz: meta.sz,
                        cn: csv.countLines()
                    ],
                    data
            )
        }
        .set {ch_reporting_cram}

    emit:
    mcool               = COOLER_ZOOMIFY.out.mcool
    ch_reporting        = ch_reporting_cram.collect()
    versions            = ch_versions.ifEmpty(null)
}

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.cram")

    "true"
}
