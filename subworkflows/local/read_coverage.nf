#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { BEDTOOLS_BAMTOBED                             } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_GENOMECOV                            } from '../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MAX          } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MIN          } from '../../modules/nf-core/bedtools/merge/main'
include { GNU_SORT                                      } from '../../modules/nf-core/gnu/sort/main'
include { MINIMAP2_ALIGN                                } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                                } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT                                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_FILTER_PRIMARY } from '../../modules/nf-core/samtools/view/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_NORMAL        } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_LOG           } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_AVGCOV           } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GRAPHOVERALLCOVERAGE                          } from '../../modules/local/graphoverallcoverage'
include { GETMINMAXPUNCHES                              } from '../../modules/local/getminmaxpunches'
include { FINDHALFCOVERAGE                              } from '../../modules/local/findhalfcoverage'
include { LONGREADCOVERAGESCALELOG                      } from '../../modules/local/longreadcoveragescalelog'
include { AVGCOV                                        } from '../../modules/local/avgcov'

workflow READ_COVERAGE {

    take:
    reference_ch        // Channel: tuple [ val(meta), file( reference_file ) ]
    dot_genome          // Channel: tuple [ val(meta), [ file( datafile ) ]   ]
    read_ch             // Channel: tuple [ val(meta), val( str )             ]  read channel (.fasta.gz)

    main:
    ch_versions                 = Channel.empty()

    //
    // LOGIC: TAKE THE READ FOLDER AS INPUT AND GENERATE THE CHANNEL OF READ FILES
    //
    ch_grabbed_reads_path       = GrabFiles( read_ch )

    ch_grabbed_reads_path
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_reads_path }

    //
    // LOGIC: PREPARE FOR MINIMAP2, USING READ_TYPE AS FILTER TO DEFINE THE MAPPING METHOD, CHECK YAML_INPUT.NF
    //
    reference_ch
        .combine( ch_reads_path )
        .combine( read_ch)
        .map { meta, ref, reads_path, read_meta, readfolder ->
            tuple(
                [   id          : meta.id,
                    single_end  : read_meta.single_end,
                    readtype    : read_meta.read_type.toString()
                ],
                reads_path,
                ref,
                true,
                false,
                false,
                read_meta.read_type.toString()
            )
        }
        .set { pre_minimap_input }

    pre_minimap_input
        .multiMap { meta, reads_path, ref, bam_output, cigar_paf, cigar_bam, reads_type ->
            read_tuple          : tuple( meta, reads_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { minimap_input }

    //
    // PROCESS: MINIMAP ALIGNMENT
    //
    MINIMAP2_ALIGN (
            minimap_input.read_tuple,
            minimap_input.ref,
            minimap_input.bool_bam_ouput,
            minimap_input.bool_cigar_paf,
            minimap_input.bool_cigar_bam
    )
    ch_versions                 = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_bams                     = MINIMAP2_ALIGN.out.bam

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
    // MODULE: MERGE ALL OUTPUT BAM
    //
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_ch,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // MODULE: SORT MAPPED BAM
    //
    SAMTOOLS_SORT (
        SAMTOOLS_MERGE.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

    //
    // MODULE: INDEXING SORTED MAPPED BAM
    //
    //SAMTOOLS_INDEX (
    //    SAMTOOLS_SORT.out.bam
    //)
    //ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    SAMTOOLS_SORT.out.bam
        .combine( reference_ch )
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
    SAMTOOLS_VIEW_FILTER_PRIMARY(
        view_input.bam_input,
        view_input.ref_input,
        []
    )
    ch_versions             = ch_versions.mix(SAMTOOLS_VIEW_FILTER_PRIMARY.out.versions)

    //
    // MODULE: BAM TO PRIMARY BED
    //
    BEDTOOLS_BAMTOBED(
        SAMTOOLS_VIEW_FILTER_PRIMARY.out.bam
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

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
    ch_versions             = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    //
    // MODULE: SORT THE PRIMARY BED FILE
    //
    GNU_SORT(
        BEDTOOLS_GENOMECOV.out.genomecov
    )
    ch_versions             = ch_versions.mix(GNU_SORT.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    GETMINMAXPUNCHES(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(GETMINMAXPUNCHES.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MAX(
        GETMINMAXPUNCHES.out.max
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_MERGE_MAX.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MIN(
        GETMINMAXPUNCHES.out.min
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_MERGE_MIN.out.versions)

    //
    // MODULE: GENERATE DEPTHGRAPH
    //
    GRAPHOVERALLCOVERAGE(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(GRAPHOVERALLCOVERAGE.out.versions)
    ch_depthgraph           = GRAPHOVERALLCOVERAGE.out.part

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
    ch_versions             = ch_versions.mix(FINDHALFCOVERAGE.out.versions)

    //
    // LOGIC: PREPARING NORMAL COVERAGE INPUT
    //
    GNU_SORT.out.sorted
        .combine( dot_genome )
        .combine(reference_ch)
        .multiMap { meta, file, meta_my_genome, my_genome, ref_meta, ref ->
            ch_coverage_bed :   tuple ([ id: ref_meta.id, single_end: true], file)
            genome_file     :   my_genome
        }
        .set { bed2bw_normal_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR NORMAL COVERAGE
    //
    BED2BW_NORMAL(
        bed2bw_normal_input.ch_coverage_bed,
        bed2bw_normal_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_NORMAL.out.versions)

    //
    // MODULE: CONVERT COVERAGE TO LOG
    //
    LONGREADCOVERAGESCALELOG(
        GNU_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix(LONGREADCOVERAGESCALELOG.out.versions)

    //
    // LOGIC: PREPARING LOG COVERAGE INPUT
    //
    LONGREADCOVERAGESCALELOG.out.bed
        .combine( dot_genome )
        .combine(reference_ch)
        .multiMap { meta, file, meta_my_genome, my_genome, ref_meta, ref ->
            ch_coverage_bed :   tuple ([ id: ref_meta.id, single_end: true], file)
            genome_file     :   my_genome
        }
        .set { bed2bw_log_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR LOG COVERAGE
    //
    BED2BW_LOG(
        bed2bw_log_input.ch_coverage_bed,
        bed2bw_log_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_LOG.out.versions)

    //
    // MODULE: CALCULATE AVERAGE COVERAGE BASED ON SCAFFOLD
    //
    AVGCOV(
        GNU_SORT.out.sorted,
        bed2bw_log_input.genome_file
    )
    ch_versions             = ch_versions.mix(AVGCOV.out.versions)

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR AVERAGE COVERAGE
    //
    BED2BW_AVGCOV(
        AVGCOV.out.avgbed,
        bed2bw_log_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_AVGCOV.out.versions)

    //
    // LOGIC: GENERATE A SUMMARY TUPLE FOR OUTPUT
    //
    ch_grabbed_reads_path
            .collect()
            .map { meta, fasta ->
                tuple( [    id: 'read',
                            sz: fasta instanceof ArrayList ? fasta.collect { it.size()} : fasta.size() ],
                            fasta
                )
            }
            .set { ch_reporting_pacbio }

    emit:
    ch_minbed               = BEDTOOLS_MERGE_MIN.out.bed
    ch_halfbed              = FINDHALFCOVERAGE.out.bed
    ch_maxbed               = BEDTOOLS_MERGE_MAX.out.bed
    ch_reporting            = ch_reporting_pacbio.collect()
    ch_covbw_nor            = BED2BW_NORMAL.out.bigwig
    ch_covbw_log            = BED2BW_LOG.out.bigwig
    ch_covbw_avg            = BED2BW_AVGCOV.out.bigwig
    versions                = ch_versions
}

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.{fa,fasta}.{gz}")

    "true"
}
