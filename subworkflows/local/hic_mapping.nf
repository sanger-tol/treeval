#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

nextflow.enable.dsl=2

// MODULE IMPORT
include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BWAMEM2_INDEX                           } from '../../modules/nf-core/bwamem2/index/main'
include { COOLER_CLOAD      } from '../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY    } from '../../modules/nf-core/cooler/zoomify/main'
include { GNU_SORT          } from '../../modules/nf-core/gnu/sort/main'
include { PRETEXTMAP as PRETEXTMAP_LOWRES  } from '../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_HIGHRES } from '../../modules/nf-core/pretextmap/main'
include { SAMTOOLS_FAIDX    } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MARKDUP  } from '../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_MERGE    } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW     } from '../../modules/nf-core/samtools/view/main'

include { GENERATE_CRAM_CSV                       } from '../../modules/local/generate_cram_csv'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT  } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { JUICER_TOOLS_PRE          } from '../../modules/local/juicer_tools_pre'
include { GET_PAIRED_CONTACT_BED    } from '../../modules/local/get_paired_contact_bed'


workflow HIC_MAPPING {
    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    reference_index     // Channel [ val(meta), path(file) ]
    dot_genome          // Channel: [val(meta), [ datafile ]]
    hic_reads_path      // Channel [ val(meta), path(directory) ]

    main:
    ch_versions         = Channel.empty()

    ch_cool_bin = Channel.of(1)

    //
    // MODULE: Indexing on reference output the folder of indexing files
    //
    BWAMEM2_INDEX (reference_tuple)
    ch_versions         = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    //
    // LOGIC: make channel of hic reads as input for GENERATE_CRAM_CSV
    //
    reference_tuple
        .combine( hic_reads_path )
        .map { meta, ref, hic_reads_path ->
                tuple([ id: meta.id, single_end: true], hic_reads_path) }
        .set { get_reads_input }

    ch_grab  = GrabFiles(get_reads_input)

    //
    // MODULE: generate a cram csv file containing the required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV ( ch_grab )
    ch_versions = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)

    //
    // LOGIC: organise all parametres into a channel for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    ch_filtering_input  = GENERATE_CRAM_CSV.out.csv
                            .splitCsv()
                            .combine (reference_tuple)
                            .combine (BWAMEM2_INDEX.out.index)
                            .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
                                  tuple([ 
                                        id: cram_id.id
                                        ], 
                                    file(cram_info[0]),
                                    cram_info[1],
                                    file(ref_dir),
                                    cram_info[2],
                                    cram_info[3],
                                    cram_info[4],
                                    cram_info[5],
                                    cram_info[6],
                                    bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1])                          
                            }

    //
    // MODULE: parallel proccessing bwa-mem2 alignment by given interval of containers from cram files
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT ( ch_filtering_input  )
    ch_versions = ch_versions.mix(CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions)

    //
    // LOGIC: PREPARING MERGE INPUT
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam
        .combine( reference_tuple )
        .combine( reference_index )
        .multiMap { bam_meta, bam, ref_meta, ref_fa, ref_idx_meta, ref_idx ->
            input_bam:  tuple(bam_meta, bam)
            reference:  ref_fa
            ref_idx:  ref_idx
        }
        .set { merge_input }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MERGE ( merge_input.input_bam, merge_input.reference, merge_input.ref_idx )
    ch_versions = ch_versions.mix ( SAMTOOLS_MERGE.out.versions.first() )

    //
    // LOGIC: PREPARING PRETEXT MAP INPUT
    //
    SAMTOOLS_MERGE.out.bam
        .combine( reference_tuple )
        .multiMap { bam_meta, bam, ref_meta, ref_fa ->
            input_bam:  tuple(bam_meta, bam)
            reference:  ref_fa
        }
        .set { pretext_input }

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR LOW RES
    //
    PRETEXTMAP_LOWRES ( pretext_input.input_bam, pretext_input.reference )
    ch_versions         = ch_versions.mix(PRETEXTMAP_LOWRES.out.versions)

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
    //
    PRETEXTMAP_HIGHRES ( pretext_input.input_bam, pretext_input.reference )
    ch_versions         = ch_versions.mix(PRETEXTMAP_HIGHRES.out.versions)

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MARKDUP ( pretext_input.input_bam, pretext_input.reference )
    ch_versions = ch_versions.mix ( SAMTOOLS_MARKDUP.out.versions.first() )

    //
    // LOGIC: PREPARING MERGED INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    SAMTOOLS_MARKDUP.out.bam
        .combine( reference_tuple )
        .combine( BWAMEM2_INDEX.out.index )
        .map { meta, file, ref_meta, ref, ref_index_meta, ref_index ->
                tuple([ id: meta.id, single_end: true], file, ref, ref_index) }
        .set { view_input }

    //
    // MODULE: GET PRIMARY BAM
    //
    SAMTOOLS_VIEW(
        view_input.map { [it[0], it[1], it[3]] },
        view_input.map { it[2] },
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: BAM TO PRIMARY BED
    //
    BEDTOOLS_BAMTOBED(SAMTOOLS_VIEW.out.bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    //
    // MODULE: SORT THE PRIMARY BED FILE
    //
    GNU_SORT(BEDTOOLS_BAMTOBED.out.bed)
    ch_versions = ch_versions.mix(GNU_SORT.out.versions)
    ch_bed = GNU_SORT.out.sorted

    //
    // MODULE: GENERATE CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED(ch_bed)
    ch_versions = ch_versions.mix(GET_PAIRED_CONTACT_BED.out.versions)

    //
    // LOGIC: PREPARE JUICER TOOLS INPUT
    //   
    GET_PAIRED_CONTACT_BED.out.bed
        .combine( dot_genome )
        .map { meta, paired_contacts, meta_my_genome, my_genome ->
                tuple([ id: meta.id, single_end: true], paired_contacts, my_genome, meta.id) }
        .set { ch_juicer_input }

    //
    // MODULE: GENERATE HIC MAP
    //
    JUICER_TOOLS_PRE(
        ch_juicer_input.map { [it[0], it[1]] },
        ch_juicer_input.map { it[2] }, 
        ch_juicer_input.map { it[3] }
    )
    ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)

    //
    // LOGIC: BIN CONTACT PAIRS
    // 
    GET_PAIRED_CONTACT_BED.out.bed
        .join(ch_bed)
        .combine(ch_cool_bin)
        .set { ch_binned_pairs }

    //
    // LOGIC: PREPARE COOLER INPUT
    //     
    ch_binned_pairs
        .combine(dot_genome)
        .map{ meta, pairs, bed, cool_bin, meta_my_genome, my_genome -> [meta, pairs, bed, cool_bin, my_genome]}
        .set { ch_cooler_input }

    //
    // MODULE: 
    //     
    COOLER_CLOAD(
        ch_cooler_input.map { [it[0], it[1], it[2], it[3]] },
        ch_cooler_input.map { it[4] }
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)
    
    //
    // LOGIC: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    //     
    COOLER_CLOAD.out.cool
        .map{ meta, cools, cool_bin -> [meta, cools]}
        .set{ch_cool}

    //
    // MODULE: ZOOM COOL TO MCOOL
    // 
    COOLER_ZOOMIFY(ch_cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    emit:
    hr_pretext      = PRETEXTMAP_HIGHRES.out.pretext
    normal_pretext  = PRETEXTMAP_LOWRES.out.pretext
    mcool           = COOLER_ZOOMIFY.out.mcool
    hic             = JUICER_TOOLS_PRE.out.hic
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
