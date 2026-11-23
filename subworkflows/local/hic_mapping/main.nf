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
    _ch_gap_file // Channel: tuple [ val(meta), path(file) ]
    _ch_coverage_file // Channel: tuple [ val(meta), path(file) ]
    _ch_telo_file // Channel: tuple [ val(meta), path(file) ]
    _ch_repeat_density_file // Channel: tuple [ val(meta), path(file) ]
    _val_workflow_setting // string: Run mode (FULL, RAPID, RAPID_TOL, etc.)
    val_create_binfile // boolean: Generate bin file using YAHS
    _val_run_juicer // boolean: Generate .hic file using Juicer
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
    mcool               = COOLER_ZOOMIFY.out.mcool
    versions            = ch_versions
}
