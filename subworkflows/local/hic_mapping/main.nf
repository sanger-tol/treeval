#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

//
// MODULE IMPORT BLOCK
//
include { COOLER_CLOAD                          } from '../../../modules/nf-core/cooler/cload/main'
include { COOLER_ZOOMIFY                        } from '../../../modules/nf-core/cooler/zoomify/main'
include { PRETEXTMAP as PRETEXTMAP_STANDRD      } from '../../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_HIGHRES      } from '../../../modules/nf-core/pretextmap/main'
include { PRETEXTMAP as PRETEXTMAP_ULTRA        } from '../../../modules/nf-core/pretextmap/main'
include { PRETEXTSNAPSHOT as SNAPSHOT_SRES      } from '../../../modules/nf-core/pretextsnapshot/main'
include { JUICERTOOLS_PRE                       } from '../../../modules/nf-core/juicertools/pre/main'
include { SUBSAMPLE_BAM                         } from '../../../modules/local/subsample/bam/main'
include { PRETEXT_GRAPH as PRETEXT_INGEST_SNDRD } from '../../../modules/local/pretext/graph/main'
include { PRETEXT_GRAPH as PRETEXT_INGEST_HIRES } from '../../../modules/local/pretext/graph/main'
include { PRETEXT_GRAPH as PRETEXT_INGEST_ULTRA } from '../../../modules/local/pretext/graph/main'
include { YAHS                                  } from '../../../modules/nf-core/yahs/main'

//
// SUBWORKFLOW IMPORT BLOCK
//
include { HIC_BAMTOBED as HIC_BAMTOBED_COOLER   } from '../hic_bamtobed/main'
include { HIC_BAMTOBED as HIC_BAMTOBED_JUICER   } from '../hic_bamtobed/main'
include { CRAM_MAP_ILLUMINA_HIC                 } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic/'

workflow HIC_MAPPING {
    take:
    ch_reference_tuple // Channel: tuple [ val(meta), path(file) ]
    ch_reference_index // Channel: tuple [ val(meta), path(file) ]
    ch_dot_genome // Channel: tuple [ val(meta), path(datafile)  ]
    ch_hic_reads_path // Channel: tuple [ val(meta), path(directory) ]
    ch_gap_file // Channel: tuple [ val(meta), path(file) ]
    ch_coverage_file // Channel: tuple [ val(meta), path(file) ]
    ch_telo_file // Channel: tuple [ val(meta), path(file) ]
    ch_repeat_density_file // Channel: tuple [ val(meta), path(file) ]
    val_workflow_setting // string: Run mode (FULL, RAPID, RAPID_TOL, etc.)
    val_create_binfile // boolean: Generate bin file using YAHS
    val_run_juicer // boolean: Generate .hic file using Juicer
    val_aligner // str: which hic aliner to use: <bwamem2/minimap2>

    main:
    ch_versions = channel.empty()

    //
    // COMMENT: 1000bp BIN SIZE INTERVALS FOR CLOAD
    //
    ch_cool_bin = channel.of(1000)

    //
    // Subworkflow: Chunked mapping of Hi-C reads to the reference with either bwamem2 or minimap2.
    // BAM has duplicates marked.
    //
    ch_cram_map_illumina_hic_input = ch_reference_tuple
        .combine(ch_hic_reads_path)
        .multiMap { meta, ref, _hic_meta, hic_cram ->
            reference: tuple(meta, ref)
            cram: tuple(meta, hic_cram)
        }

    CRAM_MAP_ILLUMINA_HIC(
        ch_cram_map_illumina_hic_input.reference,
        ch_cram_map_illumina_hic_input.cram,
        val_aligner,
        params.hic_mapping_chunk_size,
    )

    //
    // LOGIC: PREPARING PRETEXT MAP INPUT
    //
    pretext_input = CRAM_MAP_ILLUMINA_HIC.out.bam
        .combine(ch_reference_tuple)
        .combine(ch_dot_genome)
        .multiMap { bam_meta, bam, ref_meta, ref_fa, _genome_meta, genome_file ->
            input_bam: tuple([id: bam_meta.id, sz: file(bam).size()], bam)
            reference: tuple(ref_meta, ref_fa, genome_file)
        }

    //
    // LOGIC: MAKE YAHS INPUT AND VALIDATE/FIX REF/INDEX PREFIXES
    //
    ch_yahs_input = ch_reference_tuple
        .filter { val_create_binfile }
        .combine(ch_reference_index)
        .map { ref_meta, ref, _fai_meta, fai ->
            def ref_name = ref.getName()
            def expected_fai = file("${fai.parent}/${ref_name}.fai")

            if (fai.getName() == expected_fai) {
                return [ref_meta, ref, fai]
            }
            else {
                // OTHER  METHODS WERE CAUSING CHANNEL POLLUTION
                // WHERE NEW FILE NAME WOULD BE ADDED TO THE INPUT CHANNEL
                // AND CRASH ON L156
                def copy_to_dir = "${fai.parent}/renamed"
                def new_path = "${copy_to_dir}/${ref_name}.fai"

                if (!file(new_path).exists()) {
                    file(copy_to_dir).mkdirs()
                    fai.mklink(new_path)
                }

                return [ref_meta, ref, file(new_path)]
            }
        }
        .combine(CRAM_MAP_ILLUMINA_HIC.out.bam)
        .map { ref_meta, ref, fai, _bam_ref, merged_bam_path ->
            tuple(
                ref_meta,
                ref,
                fai,
                merged_bam_path,
                []) // Placeholder for AGP file input if needed in the future
        }

    //
    // MODULE: RUN YAHS TO GENERATE ALIGNMENT BIN FILE
    //
    YAHS(
        ch_yahs_input
    )

    //
    // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR LOW RES
    //
    PRETEXTMAP_STANDRD(
        pretext_input.input_bam,
        pretext_input.reference,
    )

    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //
    PRETEXT_INGEST_SNDRD(
        PRETEXTMAP_STANDRD.out.pretext,
        ch_gap_file.map { _meta, gapfile -> gapfile },
        ch_coverage_file.map { _meta, covfile -> covfile },
        ch_telo_file,
        ch_repeat_density_file.map { _meta, rdfile -> rdfile },
        params.split_telomere,
    )
    ch_versions = ch_versions.mix(PRETEXT_INGEST_SNDRD.out.versions)

    if (params.run_hires) {
        //
        // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
        //
        PRETEXTMAP_HIGHRES(
            pretext_input.input_bam,
            pretext_input.reference,
        )

        //
        // NOTICE: This could fail on LARGE hires maps due to some memory parameter in the C code
        //         of pretext graph. There is a "fixed" version in sanger /software which may need
        //         to be released in this case
        //
        // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
        //

        PRETEXT_INGEST_HIRES(
            PRETEXTMAP_HIGHRES.out.pretext,
            ch_gap_file.map { _meta, gapfile -> gapfile },
            ch_coverage_file.map { _meta, covfile -> covfile },
            ch_telo_file,
            ch_repeat_density_file.map { _meta, rdfile -> rdfile },
            params.split_telomere,
        )
        ch_versions = ch_versions.mix(PRETEXT_INGEST_HIRES.out.versions)
        hires_pretext = PRETEXT_INGEST_HIRES.out.pretext
    }
    else {
        hires_pretext = channel.empty()
    }

    if (params.run_ultra) {
        //
        // MODULE: GENERATE PRETEXT MAP FROM MAPPED BAM FOR HIGH RES
        //
        PRETEXTMAP_ULTRA(
            pretext_input.input_bam,
            pretext_input.reference,
        )

        //
        // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
        //

        PRETEXT_INGEST_ULTRA(
            PRETEXTMAP_ULTRA.out.pretext,
            ch_gap_file.map { _meta, gapfile -> gapfile },
            ch_coverage_file.map { _meta, covfile -> covfile },
            ch_telo_file,
            ch_repeat_density_file.map { _meta, rdfile -> rdfile },
            params.split_telomere,
        )
        ch_versions = ch_versions.mix(PRETEXT_INGEST_ULTRA.out.versions)
        ultra_pretext = PRETEXT_INGEST_ULTRA.out.pretext
    }
    else {
        ultra_pretext = channel.empty()
    }

    //
    // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    //
    SNAPSHOT_SRES(
        PRETEXTMAP_STANDRD.out.pretext.map { meta, pretext -> tuple(meta, pretext, []) }
    )

    //
    // LOGIC: PREPARE BAMTOBED JUICER INPUT.
    //        BRANCH TO SUBSAMPLE BAM IF LARGER THAN 50G
    //
    if (val_workflow_setting != "RAPID_TOL" && !val_run_juicer) {

        ch_merged_bam = CRAM_MAP_ILLUMINA_HIC.out.bam.branch { meta, bam ->
            def bam_sz = file(bam).size()
            tosubsample: bam_sz >= 50000000000
            return [[id: meta.id, sz: bam_sz], bam]
            unmodified: bam_sz < 50000000000
            return [[id: meta.id, sz: bam_sz], bam]
        }

        //
        // MODULE: SUBSAMPLE BAM
        //
        SUBSAMPLE_BAM(
            ch_merged_bam.tosubsample
        )
        ch_versions = ch_versions.mix(SUBSAMPLE_BAM.out.versions)

        //
        // LOGIC: COMBINE BRANCHED TO SINGLE OUTPUT
        //
        ch_subsampled_bam = SUBSAMPLE_BAM.out.subsampled_bam
        ch_subsampled_bam.mix(ch_merged_bam.unmodified)

        //
        // SUBWORKFLOW: BAM TO BED FOR JUICER - USES THE SUBSAMPLED MERGED BAM
        //
        HIC_BAMTOBED_JUICER(ch_subsampled_bam)
        ch_versions = ch_versions.mix(HIC_BAMTOBED_JUICER.out.versions)

        //
        // LOGIC: PREPARE JUICER TOOLS INPUT
        //
        ch_juicer_input = HIC_BAMTOBED_JUICER.out.paired_contacts_bed
            .combine(ch_dot_genome)
            .multiMap { meta, paired_contacts, meta_my_genome, my_genome ->
                paired: tuple([id: meta.id, single_end: true], paired_contacts)
                genome: tuple(meta_my_genome, "", my_genome)
            }


        //
        // MODULE: GENERATE HIC MAP, ONLY IS PIPELINE IS RUNNING ON MODE FULL
        //
        JUICERTOOLS_PRE(
            ch_juicer_input.paired,
            ch_juicer_input.genome,
        )
    }

    //
    // SUBWORKFLOW: BAM TO BED FOR COOLER
    //
    HIC_BAMTOBED_COOLER(CRAM_MAP_ILLUMINA_HIC.out.bam)
    ch_versions = ch_versions.mix(HIC_BAMTOBED_COOLER.out.versions)

    //
    // LOGIC: BIN CONTACT PAIRS
    //
    ch_binned_pairs = HIC_BAMTOBED_COOLER.out.paired_contacts_bed.join(HIC_BAMTOBED_COOLER.out.sorted_bed)

    //
    // MODULE: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    //
    COOLER_CLOAD(
        ch_binned_pairs,
        ch_dot_genome,
        "pairs",
        ch_cool_bin,
    )

    //
    // LOGIC: REFACTOR CHANNEL FOR ZOOMIFY
    //
    ch_cool = COOLER_CLOAD.out.cool.map { meta, cools ->
        tuple(meta, cools)
    }

    //
    // MODULE: ZOOM COOL TO MCOOL
    //
    COOLER_ZOOMIFY(ch_cool)

    emit:
    hires_pretext       = hires_pretext
    ultra_pretext       = ultra_pretext
    standardres_pretext = PRETEXT_INGEST_SNDRD.out.pretext
    standardres_png     = SNAPSHOT_SRES.out.image
    mcool               = COOLER_ZOOMIFY.out.mcool
    versions            = ch_versions
}
