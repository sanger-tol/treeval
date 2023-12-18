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
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_NORMAL    } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_LOG       } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GRAPHOVERALLCOVERAGE                      } from '../../modules/local/graphoverallcoverage'
include { GETMINMAXPUNCHES                          } from '../../modules/local/getminmaxpunches'
include { FINDHALFCOVERAGE                          } from '../../modules/local/findhalfcoverage'
include { LONGREADCOVERAGESCALELOG                  } from '../../modules/local/longreadcoveragescalelog'

workflow LONGREAD_COVERAGE {

    take:
    reference_tuple     // Channel: tuple [ val(meta), file( reference_file ) ]
    dot_genome          // Channel: tuple [ val(meta), [ file( datafile ) ]   ]
    reads_path          // Channel: tuple [ val(meta), val( str )             ]

    main:
    ch_versions             = Channel.empty()

    //
    // LOGIC: CHECK IF THE INPUT READ FILE IS PAIRED END OR SINGLE END BASED ON THE READ PLATFORM, THEN RUN MINIMAP
    //
    if ( platform.filter { it == "hifi" } || platform.filter { it == "clr" } || platform.filter { it == "ont" } ) { 
        SE_MAPPING (
            reference_tuple,
            assembly_path,
            pacbio_tuple,
            platform
        )
        ch_versions = ch_versions.mix(SE_MAPPING.out.versions)
        ch_align_bam
            .mix( SE_MAPPING.out.mapped_bam )
            .set { merged_bam }
    }
    else if ( platform.filter { it == "illumina" } ) { 

        PE_MAPPING  (
            reference_tuple,
            assembly_path,
            pacbio_tuple,
            platform
        )
        ch_versions = ch_versions.mix(PE_MAPPING.out.versions)
        ch_align_bam
            .mix( PE_MAPPING.out.mapped_bam )
            .set { merged_bam }
    }

    //
    // MODULE: SORT MAPPED BAM
    //
    SAMTOOLS_SORT (
        merged_bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT.out.versions )

    //
    // MODULE: INDEXING SORTED MAPPED BAM
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )
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
    ch_versions             = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: BAM TO PRIMARY BED
    //
    BEDTOOLS_BAMTOBED(
        SAMTOOLS_VIEW.out.bam
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
        .combine(reference_tuple)
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
        .combine(reference_tuple)
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
    // LOGIC: GENERATE A SUMMARY TUPLE FOR OUTPUT
    //
    ch_grabbed_read_paths
            .collect()
            .map { meta, fasta ->
                tuple( [    id: 'pacbio',
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
    versions                = ch_versions
}

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.fasta.gz")

    "true"
}
