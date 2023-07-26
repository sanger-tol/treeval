#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { BEDTOOLS_BAMTOBED                         } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_GENOMECOV                        } from '../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MAX      } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MIN      } from '../../modules/nf-core/bedtools/merge/main'
include { GNU_SORT                                  } from '../../modules/nf-core/gnu/sort/main'
include { MINIMAP2_INDEX                            } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SPLIT    } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN                            } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                            } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT                             } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW                             } from '../../modules/nf-core/samtools/view/main'
include { UCSC_BEDGRAPHTOBIGWIG                     } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GRAPHOVERALLCOVERAGE                      } from '../../modules/local/graphoverallcoverage'
include { GETMINMAXPUNCHES                          } from '../../modules/local/getminmaxpunches'
include { FINDHALFCOVERAGE                          } from '../../modules/local/findhalfcoverage'

// less /nfs/team135/yy5/docker_cov/run-coverage

workflow LONGREAD_COVERAGE {

    take:
    reference_tuple     // Channel: [ val(meta), file( reference_file ) ]
    dot_genome          // Channel: [ val(meta), [ file( datafile ) ]   ]
    reads_path          // Channel: [ val(meta), val( str )             ]

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: CREATES INDEX OF REFERENCE FILE
    //
    MINIMAP2_INDEX(
        reference_tuple
    )
    ch_versions = ch_versions.mix( MINIMAP2_INDEX.out.versions )

    //
    // LOGIC: PREPARE GET_READS_FROM_DIRECTORY INPUT
    //
    reference_tuple
        .combine( reads_path )
        .map { meta, ref, reads_path ->
            tuple(
                [   id          : meta.id,
                    single_end  : true  ],
                reads_path
            )
        }
        .set { get_reads_input }

    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_read_paths       = GrabFiles( get_reads_input )

    //
    // LOGIC: PACBIO READS FILES TO CHANNEL
    //
    ch_grabbed_read_paths
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_read_paths }

    //
    // LOGIC: COMBINE PACBIO READ PATHS WITH MINIMAP2_INDEX OUTPUT
    //
    MINIMAP2_INDEX.out.index
        .combine( ch_read_paths )
        .combine( reference_tuple )
        .map { meta, ref_mmi, read_path, ref_meta, reference ->
            tuple(
                [   id          : meta.id,
                    single_end  : true,
                    split_prefix: read_path.toString().split('/')[-1].split('.fasta.gz')[0]
                ],
                read_path,
                ref_mmi,
                true,
                false,
                false,
                file( reference ).size()
            )
        }
        .branch {
            large               : it[6] > 3000000000
            small               : it[6] < 3000000000
        }
        .set { mma_input }

    mma_input.large
        .multiMap { meta, read_path, ref_mmi, bam_output, cigar_paf, cigar_bam, file_size ->
            read_tuple          : tuple( meta, read_path)
            mmi_index           : ref_mmi
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { large }

    mma_input.small
        .multiMap { meta, read_path, ref_mmi, bam_output, cigar_paf, cigar_bam, file_size ->
            read_tuple          : tuple( meta, read_path)
            mmi_index           : ref_mmi
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { small }

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE <5GB PER SCAFFOLD
    //
    MINIMAP2_ALIGN (
        small.read_tuple,
        small.mmi_index,
        small.bool_bam_ouput,
        small.bool_cigar_paf,
        small.bool_cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_align_bams = MINIMAP2_ALIGN.out.bam

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE >5GB PER SCAFFOLD
    //
    MINIMAP2_ALIGN_SPLIT (
        large.read_tuple,
        large.mmi_index,
        large.bool_bam_ouput,
        large.bool_cigar_paf,
        large.bool_cigar_bam
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_SPLIT.out.versions)

    //
    // LOGIC: COLLECT OUTPUTTED BAM FILES FROM BOTH PROCESSES
    //
    ch_align_bams
        .mix( MINIMAP2_ALIGN_SPLIT.out.bam )
        .set { ch_bams }

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    ch_bams
        .map { meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [ id    : file[0].toString().split('/')[-1].split('_')[0] ], // Change sample ID
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGES THE BAM FILES IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_tuple,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // MODULE: SORT THE MERGED BAM BEFORE CONVERSION
    //
    SAMTOOLS_SORT (
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_MERGE.out.versions )

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    SAMTOOLS_SORT.out.bam
        .combine( reference_tuple )
        .multiMap { meta, bam, ref_meta, ref ->
                bam_input       :   tuple(
                                        [   id          : meta.id,
                                            sz          : bam.size(),
                                            single_end  : true  ],
                                        bam,
                                        []   // As we aren't using an index file here
                                    )
                ref_input       :   tuple(
                                        ref_meta,
                                        ref
                                    )
        }
        .set { view_input }
    //
    // MODULE: EXTRACT READS FOR PRIMARY ASSEMBLY
    //
    SAMTOOLS_VIEW(
        view_input.bam_input,
        view_input.ref_input,
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: BAM TO PRIMARY BED
    //
    BEDTOOLS_BAMTOBED(
        SAMTOOLS_VIEW.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    //
    // LOGIC: PREPARING Genome2Cov INPUT
    //
    BEDTOOLS_BAMTOBED.out.bed
        .combine( dot_genome )
        .multiMap { meta, file, my_genome_meta, my_genome ->
            input_tuple         :   tuple (
                                        [   id          :   meta.id,
                                            single_end  :   true    ],
                                        file,
                                        1
                                    )
            dot_genome          :   my_genome
            file_suffix         :   'bed'
        }
        .set { genomecov_input }

    //
    // MODULE: Genome2Cov
    //
    BEDTOOLS_GENOMECOV(
        genomecov_input.input_tuple,
        genomecov_input.dot_genome,
        genomecov_input.file_suffix
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    //
    // MODULE: SORT THE PRIMARY BED FILE
    //
    GNU_SORT(
        BEDTOOLS_GENOMECOV.out.genomecov
    )
    ch_versions = ch_versions.mix(GNU_SORT.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    GETMINMAXPUNCHES(
        GNU_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(GETMINMAXPUNCHES.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MAX(
        GETMINMAXPUNCHES.out.max
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_MAX.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MIN(
        GETMINMAXPUNCHES.out.min
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_MIN.out.versions)

    //
    // MODULE: GENERATE DEPTHGRAPH
    //
    GRAPHOVERALLCOVERAGE(
        GNU_SORT.out.sorted
    )
    ch_versions = ch_versions.mix(GRAPHOVERALLCOVERAGE.out.versions)
    ch_depthgraph = GRAPHOVERALLCOVERAGE.out.part

    //
    // LOGIC: PREPARING FINDHALFCOVERAGE INPUT
    //
    GNU_SORT.out.sorted
        .combine( GRAPHOVERALLCOVERAGE.out.part )
        .combine( dot_genome )
        .multiMap { meta, file, meta_depthgraph, depthgraph, meta_my_genome, my_genome ->
            halfcov_bed     :       tuple( [ id : meta.id, single_end : true  ], file )
            genome_file     :       my_genome
            depthgraph_file :       depthgraph
        }
        .set { halfcov_input }

    //
    // MODULE: FIND REGIONS OF HALF COVERAGE
    //
    FINDHALFCOVERAGE(
        halfcov_input.halfcov_bed,
        halfcov_input.genome_file,
        halfcov_input.depthgraph_file
    )
    ch_versions = ch_versions.mix(FINDHALFCOVERAGE.out.versions)

    //
    // LOGIC: PREPARING COVERAGE INPUT
    //
    GNU_SORT.out.sorted
        .combine( dot_genome )
        .multiMap { meta, file, meta_my_genome, my_genome ->
            ch_coverage_bed :   tuple ([ id: meta.id, single_end: true], file)
            genome_file     :   my_genome
        }
        .set { bed2bw_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG
    //
    UCSC_BEDGRAPHTOBIGWIG(
        bed2bw_input.ch_coverage_bed,
        bed2bw_input.genome_file
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    //
    //
    //
    ch_read_paths
        .map { meta, read_path ->
            tuple(
                [   id          : meta.id,
                    sz          : file( read_path ).size()
                ],
                read_path,
            )
        }
        .set { ch_for_PacBio }

    emit:
    ch_minbed       = BEDTOOLS_MERGE_MIN.out.bed
    ch_halfbed      = FINDHALFCOVERAGE.out.bed
    ch_maxbed       = BEDTOOLS_MERGE_MAX.out.bed
    ch_bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    ch_reporting    = ch_for_PacBio.collect()
    versions        = ch_versions
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.fasta.gz")

    "true"
}
