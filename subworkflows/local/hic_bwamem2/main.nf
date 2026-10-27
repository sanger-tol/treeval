#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

//
// MODULE IMPORT BLOCK
//
include { BWAMEM2_INDEX                                   } from '../../../modules/nf-core/bwamem2/index/main'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT          } from '../../../modules/local/cram/filter_align_bwamem2_fixmate_sort/main'
include { SAMTOOLS_MERGE                                  } from '../../../modules/nf-core/samtools/merge/main'

workflow HIC_BWAMEM2 {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path( file )      ]
    csv_ch
    reference_index

    main:
    ch_versions         = channel.empty()
    mappedbam_ch        = channel.empty()

    BWAMEM2_INDEX (
        reference_tuple
        )

    csv_ch
        .splitCsv()
        .combine ( reference_tuple )
        .combine ( BWAMEM2_INDEX.out.index )
        .map{ cram_id, cram_info, _ref_id, ref_dir, _bwa_id, bwa_path ->
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
                bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1],
                ref_dir
            )
    }
    .set { ch_filtering_input }

    //
    // MODULE: map hic reads by 10,000 container per time using bwamem2
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT (
        ch_filtering_input

    )
    ch_versions         = ch_versions.mix( CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions )
    mappedbam_ch        = CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam

    //
    // LOGIC: PREPARING BAMS FOR MERGE
    //
    mappedbam_ch
        .map{ _meta, file ->
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
    // LOGIC: PREPARING REFERENCE FOR MERGE
    //
    reference_tuple
        .combine ( reference_index )
        .map{ _ref_meta, ref, _ref_index_meta, ref_index ->
            tuple(
                [id: _ref_meta.id],
                ref,
                ref_index,
                [])
        }
        .set { reference_for_merge }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MERGE (
        collected_files_for_merge,
        reference_for_merge
    )

    emit:
    mergedbam           = SAMTOOLS_MERGE.out.bam
    versions            = ch_versions
}
