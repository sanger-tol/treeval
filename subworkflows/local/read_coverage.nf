#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { BEDTOOLS_BAMTOBED                             } from '../../modules/nf-core/bedtools/bamtobed/main'
include { BEDTOOLS_GENOMECOV                            } from '../../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MAX          } from '../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MERGE as BEDTOOLS_MERGE_MIN          } from '../../modules/nf-core/bedtools/merge/main'
include { GNU_SORT as GNU_SORT_BED                      } from '../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as GNU_SORT_COVBED                   } from '../../modules/nf-core/gnu/sort/main'
include { CAT_CAT                                       } from '../../modules/nf-core/cat/cat/main'
include { MINIMAP2_ALIGN                                } from '../../modules/nf-core/minimap2/align/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_NORMAL        } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as BED2BW_AVGCOV        } from '../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { GRAPHOVERALLCOVERAGE                          } from '../../modules/local/graphoverallcoverage'
include { GETMINMAXPUNCHES                              } from '../../modules/local/getminmaxpunches'
include { FINDHALFCOVERAGE                              } from '../../modules/local/findhalfcoverage'
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
    read_ch
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
        .combine( read_ch )
        .map { meta, ref, reads_path, read_meta, readfolder ->
            tuple(
                [   id          : meta.id,
                    single_end  : read_meta.single_end,
                    readtype    : read_meta.read_type.toString()
                ],
                reads_path,
                ref,
                false,
                false,
                false,
                true,
                read_meta.read_type.toString()
            )
        }
        .set { pre_minimap_input }

    pre_minimap_input
        .multiMap { meta, reads_path, ref, bam_output, cigar_paf, cigar_bam, bed_output, reads_type ->
            read_tuple          : tuple( meta, reads_path)
            ref                 : tuple( meta, ref)
            bool_bam_ouput      : bam_output
            val_bam_index       : "bai"
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
            bool_bed_output     : bed_output
        }
        .set { minimap_input }

    //
    // PROCESS: MINIMAP ALIGNMENT
    //
    MINIMAP2_ALIGN (
            minimap_input.read_tuple,
            minimap_input.ref,
            minimap_input.bool_bam_ouput,
            minimap_input.val_bam_index,
            minimap_input.bool_cigar_paf,
            minimap_input.bool_cigar_bam,
            minimap_input.bool_bed_output
    )
    ch_versions                 = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_beds                     = MINIMAP2_ALIGN.out.bed

    ch_beds
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
    // MODULE: MERGE ALL OUTPUT BEDS
    //
    CAT_CAT(
        collected_files_for_merge
    )
    ch_versions             = ch_versions.mix( CAT_CAT.out.versions )

    //
    // MODULE: SORT THE MERGED BED FILE INTO CHROMOSOME-LOCATION ORDER
    //
    GNU_SORT_BED(
        CAT_CAT.out.file_out
    )
    ch_versions             = ch_versions.mix(GNU_SORT_BED.out.versions)
    ch_sorted_bed           = GNU_SORT_BED.out.sorted

    //
    // LOGIC: PREPARING Genome2Cov INPUT
    //
    ch_sorted_bed
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
        genomecov_input.file_suffix,
        false
    )
    ch_versions             = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    //
    // LOGIC: BED2BIGWIG TAKES SORTED COVERAGE BED FILE
    //
    GNU_SORT_COVBED(
        BEDTOOLS_GENOMECOV.out.genomecov
    )
    ch_versions             = ch_versions.mix(GNU_SORT_COVBED.out.versions)
    ch_sorted_covbed        = GNU_SORT_COVBED.out.sorted

    //
    // MODULE: get_minmax_punches
    //
    GETMINMAXPUNCHES(
        ch_sorted_covbed
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
        ch_sorted_covbed
    )
    ch_versions             = ch_versions.mix(GRAPHOVERALLCOVERAGE.out.versions)
    ch_depthgraph           = GRAPHOVERALLCOVERAGE.out.part

    //
    // LOGIC: PREPARING FINDHALFCOVERAGE INPUT
    //
    ch_sorted_covbed
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
    ch_sorted_covbed
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
    // MODULE: CALCULATE AVERAGE COVERAGE BASED ON SCAFFOLD
    //
    AVGCOV(
        ch_sorted_covbed,
        bed2bw_normal_input.genome_file
    )
    ch_versions             = ch_versions.mix(AVGCOV.out.versions)

    //
    // MODULE: CONVERT BEDGRAPH TO BIGWIG FOR AVERAGE COVERAGE
    //
    BED2BW_AVGCOV(
        AVGCOV.out.avgbed,
        bed2bw_normal_input.genome_file
    )
    ch_versions             = ch_versions.mix(BED2BW_AVGCOV.out.versions)

    //
    // LOGIC: GENERATE A SUMMARY TUPLE FOR OUTPUT
    //
    read_ch
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
    ch_covbw_avg            = BED2BW_AVGCOV.out.bigwig
    versions                = ch_versions
}
