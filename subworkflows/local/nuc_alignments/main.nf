#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { MINIMAP2_ALIGN        } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE        } from '../../../modules/nf-core/samtools/merge/main'
include { BEDTOOLS_SORT         } from '../../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_BAMTOBED     } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { UCSC_BEDTOBIGBED      } from '../../../modules/nf-core/ucsc/bedtobigbed/main'
include { PAFTOOLS_SAM2PAF      } from '../../../modules/nf-core/paftools/sam2paf/main'
include { PAF2BED               } from '../../../modules/local/paf/to_bed/main'

//
// SUBWORKFLOW IMPORTS
//
include { PUNCHLIST             } from '../punchlist/main'

workflow NUC_ALIGNMENTS {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path(file) ]
    reference_index     // Channel: tuple [ val(meta), path(file) ]
    nuc_files           // Channel: tuple [ val(meta), path(file) ]
    dot_genome          // Channel: tuple [ val(meta), path(file) ]
    intron_size         // Channel: val(50k)

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: COLLECTION FROM GENE_ALIGNMENT IS A LIST OF ALL META AND ALL FILES
    //        BELOW CONVERTS INTO TUPLE FORMAT AND ADDS BOOLEANS FOR MINIMAP2_ALIGN
    //
    nuc_files
        .flatten()
        .buffer( size: 2 )
        .combine ( reference_tuple )
        .combine( intron_size )
        .map { meta, nuc_file, ref_meta, ref, intron ->
            tuple( [id:             meta.id,
                    type:           meta.type,
                    org:            meta.org,
                    intron_size:    intron,
                    split_prefix:   nuc_file.toString().split('/')[-1].split('.fasta')[0],
                    single_end:     true
                    ],
                    nuc_file,
                    ref,
                    true,
                    false,
                    false,
                    false
            )
        }
        .multiMap {meta, nuc_file, reference, bool_1, bool_2, bool_3, bool_4 ->
            nuc             : tuple(meta, nuc_file)
            ref             : tuple(meta, reference)
            bool_bam_output : bool_1
            val_bam_output  : "bai"
            bool_cigar_paf  : bool_2
            bool_cigar_bam  : bool_3
            bool_bedfile    : bool_4
        }
        .set { formatted_input }

    //
    // MODULE: ALIGNS REFERENCE FAIDX TO THE GENE_ALIGNMENT QUERY FILE FROM NUC_FILES
    //         EMITS ALIGNED BAM FILE
    //
    MINIMAP2_ALIGN (
        formatted_input.nuc,
        formatted_input.ref,
        formatted_input.bool_bam_output,
        formatted_input.val_bam_output,
        formatted_input.bool_cigar_paf,
        formatted_input.bool_cigar_bam,
        formatted_input.bool_bedfile
    )
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // LOGIC: CONVERTS THE MINIMAP OUTPUT TUPLE INTO A GROUPED TUPLE PER INPUT QUERY ORGANISM
    //        AND DATA TYPE (RNA, CDS, DNA).
    //
    MINIMAP2_ALIGN.out.bam
        .map { meta, file ->
            tuple(
                [   id: meta.org,
                    type: meta.type ],
                file) }
        .groupTuple( by: [0] )  // group by meta list
        .set { merge_input }

    //
    // MODULE: MERGES THE BAM FILES FOUND IN THE GROUPED TUPLE IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
    SAMTOOLS_MERGE (
        merge_input,
        reference_tuple,
        reference_index
    )
    ch_versions     = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // SUBWORKFLOW: GENERATES A PUNCHLIST FROM MERGED BAM FILE
    //
    PUNCHLIST (
        reference_tuple,
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions     = ch_versions.mix(PUNCHLIST.out.versions)

    //
    // MODULE: CONVERTS THE ABOVE MERGED BAM INTO BED FORMAT
    //
    BEDTOOLS_BAMTOBED ( SAMTOOLS_MERGE.out.bam )
    ch_versions     = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    // TODO: try filtering out here too

    //
    // LOGIC: ADDING LINE COUNT TO THE FILE FOR BETTER RESOURCE USAGE
    //
    BEDTOOLS_BAMTOBED.out.bed
        .map { meta, file ->
            tuple ( [   id:     meta.id,
                        type:   meta.type,
                        lines:  file.countLines()
                    ],
                    file
            )
        }
        .set { bedtools_input }

    //
    // MODULE: SORTS THE ABOVE BED FILE
    //
    BEDTOOLS_SORT (
        bedtools_input,
        []
    )
    ch_versions     = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    //
    // LOGIC: COMBINES GENOME_FILE CHANNEL AND ABOVE OUTPUT, SPLITS INTO TWO CHANNELS
    //        ALSO FILTERS OUT EMPTY MERGED.BED BASED ON WHETHER FILE IS >141 BYTES
    //
    BEDTOOLS_SORT.out.sorted
        .map { meta, file ->
                tuple( [    id:         meta.id,
                            type:       meta.type,
                            file_size:  file.size()
                        ],
                        file ) }
        .filter { it[0].file_size >= 141 } // Take the first item in input (meta) and check if size is more than a symlink
        .combine( dot_genome )
        .multiMap { meta, ref, genome_meta, genome ->
            bed_file:   tuple( [    id:         meta.id,
                                    type:       meta.type,
                                ],
                                ref )
            dot_genome: genome
        }
        .set { ucsc_input }

    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGBED FILE
    //
    UCSC_BEDTOBIGBED (
        ucsc_input.bed_file,
        ucsc_input.dot_genome,
        []
    )
    ch_versions     = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )

    emit:
    nuc_alignment   = UCSC_BEDTOBIGBED.out.bigbed.collect()
    punchlist       = PUNCHLIST.out.punchlist.collect()
    versions        = ch_versions
}
