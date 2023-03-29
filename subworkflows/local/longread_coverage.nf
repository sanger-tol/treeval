#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MAX } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MIN } from '../../modules/nf-core/bedtools/merge/main'
include { MINIMAP2_INDEX } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SPLIT } from '../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'

include { GET_READS_FROM_DIRECTORY } from '../../modules/local/get_reads_from_directory'
include { LONGREADCOVERAGE_GRAPHOVERALL } from '../../modules/local/longreadcoverage_graphoverall'
include { LONGREADCOVERAGE_GETMINMAXPUNCHES } from '../../modules/local/longreadcoverage_getminmaxpunches'
include { LONGREADCOVERAGE_FINDHALFCOVERAGE } from '../../modules/local/longreadcoverage_findhalfcoverage'
include { SORT_INTERSECT } from '../../modules/local/sort_intersect'

// less /nfs/team135/yy5/docker_cov/run-coverage

workflow LONGREAD_COVERAGE {

    take:
    reference_tuple     // Channel: [id: sample_id], reference_file
    dot_genome          // Channel: [val(meta), [ datafile ]]
    reads_path          // Channel: val( str )
    size_class          // Channel: val( str )

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: CREATES INDEX OF REFERENCE FILE
    //
    MINIMAP2_INDEX(reference_tuple)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    ch_ref_index = MINIMAP2_INDEX.out.index

    //
    // LOGIC: PREPARE GET_READS_FROM_DIRECTORY INPUT 
    //
    reference_tuple
        .combine( reads_path )
        .map { meta, ref, reads_path ->
                tuple([ id: meta.id, single_end: true], reads_path) }
        .set { get_reads_input }

    //
    // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
    //
    GET_READS_FROM_DIRECTORY( 
        get_reads_input.map { [it[0], it[1]] }           
    )

    //
    // LOGIC: PACBIO READS FILES TO CHANNEL
    //
    GET_READS_FROM_DIRECTORY.out.files
           .map { meta, files ->
            tuple(files)
            }
        .flatten()
        .set { ch_read_paths }

    //
    // LOGIC: COMBINE PACBIO READ PATHS WITH MINIMAP2_INDEX OUTPUT
    //
    ch_ref_index
        .combine(ch_read_paths)
        .combine(size_class)
        .map { meta, ref_mmi, read_path, size_class ->
            tuple([ id: meta.id,
                    single_end: true,
                    split_prefix: read_path.toString().split('/')[-1].split('.fasta.gz')[0]
                ],
                read_path, ref_mmi, true, false, false, size_class)
            }
        .branch {
            large: it[6] == 'L'
            small: it[6] == 'S'
        }
        .set { mma_input }

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE <5GB PER SCAFFOLD
    //   
    MINIMAP2_ALIGN (
        mma_input.small.map { [it[0], it[1]] },
        mma_input.small.map { it[2] },
        mma_input.small.map { it[3] },
        mma_input.small.map { it[4] },
        mma_input.small.map { it[5] }
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_align_bams = MINIMAP2_ALIGN.out.bam

    //
    // MODULE: ALIGN READS TO REFERENCE WHEN REFERENCE >5GB PER SCAFFOLD
    //
    MINIMAP2_ALIGN_SPLIT (
        mma_input.large.map { [it[0], it[1]] },
        mma_input.large.map { it[2] },
        mma_input.large.map { it[3] },
        mma_input.large.map { it[4] },
        mma_input.large.map { it[5] }
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_SPLIT.out.versions)
    ch_split_bams = MINIMAP2_ALIGN_SPLIT.out.bam

    //
    // LOGIC: COLLECT OUTPUTTED BAM FILES FROM BOTH PROCESSES
    //        
    ch_align_bams
        .mix(ch_split_bams)
        .set { ch_bams }

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    ch_bams
        .combine( reference_tuple )
        .combine( ch_ref_index )
        .map { meta, file, ref_meta, ref, ref_index_meta, ref_index ->
                tuple([ id: meta.id, single_end: true], file, ref, ref_index) }
        .set { merge_input }

    //
    // MODULE: MERGES THE BAM FILES IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM

    SAMTOOLS_MERGE(
        merge_input.map { [it[0], it[1]] },
        merge_input.map { it[2] }, 
        merge_input.map { it[3] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    ch_merged_bam = SAMTOOLS_MERGE.out.bam

    //
    // LOGIC: PREPARING MERGE INPUT WITH REFERENCE GENOME AND REFERENCE INDEX
    //
    ch_merged_bam
        .combine( reference_tuple )
        .combine( ch_ref_index )
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

    // //
    // // MODULE: SORTS THE PRIMARY BED FILE
    // //
    // SORT_INTERSECT ( BEDTOOLS_BAMTOBED.out.bed )

    //
    // LOGIC: PREPARING Genome2Cov INPUT
    //
    BEDTOOLS_BAMTOBED.out.bed
        .combine(dot_genome)
        .map { meta, file, my_genome_meta, my_genome -> 
            tuple([ id: meta.id, single_end: true], file, 1, my_genome, 'bed')
        }
        .set { genomecov_input }

    //
    // MODULE: Genome2Cov
    // 
    BEDTOOLS_GENOMECOV(
        genomecov_input.map { [it[0], it[1], it[2]] },
        genomecov_input.map { it[3] },
        genomecov_input.map { it[4] }
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)
    ch_coverage_unsorted_bed = BEDTOOLS_GENOMECOV.out.genomecov

    //
    // MODULE: SORTS THE PRIMARY BED FILE
    //
    SORT_INTERSECT ( ch_coverage_unsorted_bed )
    ch_coverage_bed = SORT_INTERSECT.out.sorted

    //
    // MODULE: get_minmax_punches
    //
    LONGREADCOVERAGE_GETMINMAXPUNCHES(
        ch_coverage_bed
    )
    ch_versions = ch_versions.mix(LONGREADCOVERAGE_GETMINMAXPUNCHES.out.versions)

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MAX(
        LONGREADCOVERAGE_GETMINMAXPUNCHES.out.max
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_MAX.out.versions)
    ch_maxbed = BEDTOOLS_MERGE_MAX.out.bed

    //
    // MODULE: get_minmax_punches
    //
    BEDTOOLS_MERGE_MIN(
        LONGREADCOVERAGE_GETMINMAXPUNCHES.out.min
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE_MIN.out.versions)
    ch_minbed = BEDTOOLS_MERGE_MIN.out.bed

    //
    // MODULE: GENERATE DEPTHGRAPH
    //
    LONGREADCOVERAGE_GRAPHOVERALL(
        ch_coverage_bed
    )
    ch_versions = ch_versions.mix(LONGREADCOVERAGE_GRAPHOVERALL.out.versions)
    ch_depthgraph = LONGREADCOVERAGE_GRAPHOVERALL.out.part

    //
    // LOGIC: PREPARING FINDHALFCOVERAGE INPUT
    //
    ch_coverage_bed
        .combine( ch_depthgraph )
        .combine( dot_genome )
        .map { meta, file, meta_depthgraph, depthgraph, meta_my_genome, my_genome -> 
            tuple([ id: meta.id, single_end: true], file, my_genome, depthgraph)
        }
        .set { findhalfcov_input }

    //
    // MODULE: findHalfcoverage
    //
    LONGREADCOVERAGE_FINDHALFCOVERAGE(
        findhalfcov_input.map { [it[0], it[1]] },
        findhalfcov_input.map { it[2] },
        findhalfcov_input.map { it[3] }
    )
    ch_versions = ch_versions.mix(LONGREADCOVERAGE_FINDHALFCOVERAGE.out.versions)
    ch_halfbed = LONGREADCOVERAGE_FINDHALFCOVERAGE.out.bed

    //
    // LOGIC: PREPARING FINDHALFCOVERAGE INPUT
    //
    ch_coverage_bed
        .combine( dot_genome )
        .map { meta, file, meta_my_genome, my_genome -> 
            tuple([ id: meta.id, single_end: true], file, my_genome)
        }
        .set { bed2bw_input }

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG
    //
    UCSC_BEDGRAPHTOBIGWIG(
        bed2bw_input.map { [it[0], it[1]] },
        bed2bw_input.map { it[2] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG.out.bigwig

    emit:
    ch_minbed
    ch_halfbed
    ch_maxbed
    ch_bigwig
    versions = ch_versions
}