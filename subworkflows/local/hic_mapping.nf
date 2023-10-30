#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

//
// MODULE IMPORT BLOCK
//
include { BWAMEM2_INDEX                             } from '../../modules/nf-core/bwamem2/index/main'
include { COOLER_CLOAD                              } from '../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY                            } from '../../modules/nf-core/cooler/zoomify/main'
include { PRETEXTMAP as PRETEXTMAP_STANDRD          } from '../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_HIGHRES          } from '../../modules/nf-core/pretextmap/main'
include { PRETEXTSNAPSHOT as SNAPSHOT_SRES          } from '../../modules/nf-core/pretextsnapshot/main'
include { PRETEXTSNAPSHOT as SNAPSHOT_HRES          } from '../../modules/nf-core/pretextsnapshot/main'
include { SAMTOOLS_MARKDUP                          } from '../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_MERGE                            } from '../../modules/nf-core/samtools/merge/main'
include { BAMTOBED_SORT                             } from '../../modules/local/bamtobed_sort.nf'
include { GENERATE_CRAM_CSV                         } from '../../modules/local/generate_cram_csv'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT    } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { JUICER_TOOLS_PRE                          } from '../../modules/local/juicer_tools_pre'
include { GET_PAIRED_CONTACT_BED                    } from '../../modules/local/get_paired_contact_bed'
include { PRETEXT_INGESTION as PRETEXT_INGEST_SNDRD } from '../../subworkflows/local/pretext_ingestion'
include { PRETEXT_INGESTION as PRETEXT_INGEST_HIRES } from '../../subworkflows/local/pretext_ingestion'


workflow HIC_MAPPING {
    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    reference_index     // Channel [ val(meta), path(file) ]
    dot_genome          // Channel [ val(meta), [ datafile ]]
    hic_reads_path      // Channel [ val(meta), path(directory) ]
    assembly_id         // Channel val( id )
    gap_file
    coverage_file
    logcoverage_file
    telo_file
    repeat_density_file
    workflow_setting    // val( {RAPID | FULL } )

    main:
    ch_versions         = Channel.empty()

    // COMMENT: 1000bp BIN SIZE INTERVALS FOR CLOAD
    ch_cool_bin         = Channel.of( 1000 )

    //
    // MODULE: Indexing on reference output the folder of indexing files
    //
    BWAMEM2_INDEX (
        reference_tuple
    )
    ch_versions         = ch_versions.mix( BWAMEM2_INDEX.out.versions )

    //
    // LOGIC: make channel of hic reads as input for GENERATE_CRAM_CSV
    //
    reference_tuple
        .combine( hic_reads_path )
        .map { meta, ref, hic_reads_path ->
                tuple(
                    [ id: meta.id, single_end: true],
                    hic_reads_path
                )
        }
        .set { get_reads_input }

    //
    // MODULE: generate a cram csv file containing the required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV (
        get_reads_input
    )
    ch_versions         = ch_versions.mix( GENERATE_CRAM_CSV.out.versions )

    //
    // LOGIC: organise all parametres into a channel for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV.out.csv
        .splitCsv()
        .combine (reference_tuple)
        .combine (BWAMEM2_INDEX.out.index)
        .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
                tuple([
                        id: cram_id.id
                        ],
                    file(cram_info[0]),
                    cram_info[1],
                    cram_info[2],
                    cram_info[3],
                    cram_info[4],
                    cram_info[5],
                    cram_info[6],
                    bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1]
                )
        }
        .set { ch_filtering_input }

    //
    // MODULE: parallel proccessing bwa-mem2 alignment by given interval of containers from cram files
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT (
        ch_filtering_input
    )
    ch_versions         = ch_versions.mix( CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions )

    //
    // LOGIC: PREPARING BAMS FOR MERGE
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam
        .map{ meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [
                id: file[0].toString().split('/')[-1].split('_')[0] + '_' + file[0].toString().split('/')[-1].split('_')[1]
                ],
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MERGE (
        collected_files_for_merge,
        reference_tuple,
        reference_index
    )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MERGE.out.versions.first() )

    //
    // LOGIC: PREPARING PRETEXT MAP INPUT
    //
    SAMTOOLS_MERGE.out.bam
        .combine( reference_tuple )
        .multiMap { bam_meta, bam, ref_meta, ref_fa ->
            input_bam:  tuple( [    id: bam_meta.id,
                                    sz: file( bam ).size() ],
                                bam
                        )
            reference:  ref_fa
        }
        .set { pretext_input }

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR LOW RES
    //
    PRETEXTMAP_STANDRD (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix( PRETEXTMAP_STANDRD.out.versions )

    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //
    PRETEXT_INGEST_SNDRD (
        PRETEXTMAP_STANDRD.out.pretext,
        gap_file,
        coverage_file,
        logcoverage_file,
        telo_file,
        repeat_density_file
    )
    ch_versions         = ch_versions.mix( PRETEXT_INGEST_SNDRD.out.versions )

    //
    // LOGIC: HIRES IS TOO INTENSIVE FOR RUNNING IN GITHUB CI SO THIS STOPS IT RUNNING
    //
    if ( params.config_profile_name ) {
        config_profile_name = params.config_profile_name
    } else {
        config_profile_name = 'Local'
    }

    if ( !config_profile_name.contains('GitHub') ) {
        //
        // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
        //
        PRETEXTMAP_HIGHRES (
            pretext_input.input_bam,
            pretext_input.reference
        )
        ch_versions         = ch_versions.mix( PRETEXTMAP_HIGHRES.out.versions )

        PRETEXT_INGEST_HIRES (
            PRETEXTMAP_HIGHRES.out.pretext,
            gap_file,
            coverage_file,
            logcoverage_file,
            telo_file,
            repeat_density_file
        )
        ch_versions         = ch_versions.mix( PRETEXT_INGEST_HIRES.out.versions )
    }

    //
    // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    //
    SNAPSHOT_SRES (
        PRETEXTMAP_STANDRD.out.pretext
    )
    ch_versions         = ch_versions.mix ( SNAPSHOT_SRES.out.versions )

    // NOTE: CURRENTLY UNDER INVESTIGATION
    //
    // MODULE: GENERATE PNG FROM HIGHRES PRETEXT
    //
    // SNAPSHOT_HRES ( PRETEXTMAP_HIGHRES.out.pretext )
    // ch_versions         = ch_versions.mix ( SNAPSHOT_HRES.out.versions )

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MARKDUP (
        pretext_input.input_bam,
        pretext_input.reference
    )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MARKDUP.out.versions )

    //
    // MODULE: SAMTOOLS FILTER OUT DUPLICATE READS | BAMTOBED | SORT BED FILE
    //
    BAMTOBED_SORT(
        SAMTOOLS_MARKDUP.out.bam
    )
    ch_versions         = ch_versions.mix( BAMTOBED_SORT.out.versions )

    //
    // MODULE: GENERATE CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED( BAMTOBED_SORT.out.sorted_bed )
    ch_versions         = ch_versions.mix( GET_PAIRED_CONTACT_BED.out.versions )

    //
    // LOGIC: SECTION ONLY NEEDED FOR TREEVAL VISUALISATION, NOT RAPID ANALYSIS
    //
    if (workflow_setting == 'FULL' && !config_profile_name.contains('GitHub')) {
        //
        // LOGIC: PREPARE JUICER TOOLS INPUT
        //
        GET_PAIRED_CONTACT_BED.out.bed
            .combine( dot_genome )
            .multiMap {  meta, paired_contacts, meta_my_genome, my_genome ->
                paired      :   tuple([ id: meta.id, single_end: true], paired_contacts )
                genome      :   my_genome
                id          :   meta.id
            }
            .set { ch_juicer_input }

        //
        // MODULE: GENERATE HIC MAP, ONLY IS PIPELINE IS RUNNING ON ENTRY FULL
        //

        JUICER_TOOLS_PRE(
            ch_juicer_input.paired,
            ch_juicer_input.genome,
            ch_juicer_input.id
        )
        ch_versions         = ch_versions.mix( JUICER_TOOLS_PRE.out.versions )
    }

    //
    // LOGIC: BIN CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED.out.bed
        .join( BAMTOBED_SORT.out.sorted_bed )
        .combine( ch_cool_bin )
        .set { ch_binned_pairs }

    //
    // LOGIC: PREPARE COOLER INPUT
    //
    ch_binned_pairs
        .combine(dot_genome)
        .multiMap { meta, pairs, bed, cool_bin, meta_my_genome, my_genome ->
            cooler_in   : tuple ( meta, pairs, bed, cool_bin )
            genome_file : my_genome
        }
        .set { ch_cooler }

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
        .map{ meta, cools, cool_bin ->
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

    ch_cram_files = GrabFiles( get_reads_input )

    ch_cram_files
        .collect()
        .map { meta, cram ->
            tuple( [    id: 'cram',
                        sz: cram instanceof ArrayList ? cram.collect { it.size()} : cram.size() ],
                    cram
            )
        }
        .set { ch_reporting_cram }

    emit:
    mcool               = COOLER_ZOOMIFY.out.mcool
    ch_reporting        = ch_reporting_cram.collect()
    versions            = ch_versions.ifEmpty(null)
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.cram")

    "true"
}
