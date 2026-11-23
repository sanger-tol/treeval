// SANGER-TOL SUBWORKFLOWS
include { PRETEXT_ACCESSORY_FILES                           } from '../../../subworkflows/sanger-tol/pretext_accessory_files/main'
include { CRAM_MAP_ILLUMINA_HIC                             } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic/main'
include { PAIRS_CREATE_CONTACT_MAPS as CREATE_MAPS_STDRD    } from '../../../subworkflows/sanger-tol/pairs_create_contact_maps/main'
include { PAIRS_CREATE_CONTACT_MAPS as CREATE_MAPS_HIRES    } from '../../../subworkflows/sanger-tol/pairs_create_contact_maps/main'
include { PAIRS_CREATE_CONTACT_MAPS as CREATE_MAPS_ULTRA    } from '../../../subworkflows/sanger-tol/pairs_create_contact_maps/main'

// SANGER-TOL MODULES
include { PRETEXTANNOTATE                                   } from '../../../modules/sanger-tol/pretextannotate/main'

// LOCAL MODULES
include { PRETEXT_GRAPH as PRETEXT_INGEST_SNDRD             } from '../../../modules/local/pretext/graph/main'
include { PRETEXT_GRAPH as PRETEXT_INGEST_HIRES             } from '../../../modules/local/pretext/graph/main'
include { PRETEXT_GRAPH as PRETEXT_INGEST_ULTRA             } from '../../../modules/local/pretext/graph/main'

// NF-CORE MODULES
include { SAMTOOLS_FLAGSTAT                                 } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_INDEX                                    } from '../../../modules/nf-core/samtools/index/main'

workflow CURATIONPRETEXT {
    take:
    ch_reference_tuple  // Channel [ val(meta), path(file)   ]
    ch_reference_sizes  // Channel [ val(meta), path(file)   ]
    ch_longread_reads   // Channel [ val(meta), [ path(file) ] ]
    ch_teloseq          // Channel [ val(meta), val(telomere_sequence) ]
    ch_hic_reads        // Channel [ val(meta), path(file) ]
    val_selected_aligner// val(bool)
    pre_mapped_bam      // Channel [ val(meta), path(file) ]
    val_cram_chunk_size // val(int)
    ch_snapshot_order   // channel( path(file) )
    val_snapshot_generation
    val_snapshot_annotate // val(bool)
    val_mapping_statistics // val(bool)
    val_juicer          // val(bool)
    map_params          // [ val(name): val(bool) ]

    main:
    //
    // LOGIC: IN SOME CASES THE USER MAY NOT NEED ALL OR A SELECT GROUP OF
    //          ACCESSORY FILES SO WE HAVE AN OPTION TO TURN THEM OFF
    //

    log.info "ACCESSORY TRACK OPTIONS: $map_params"

    if (!map_params["generate_tracks"]) {
        gaps_file       = channel.empty()
        cove_file       = channel.empty()
        telo_file       = channel.empty()
        rept_file       = channel.empty()
        _busc_file       = channel.empty()
        _pebb_file       = channel.empty()

    } else {
        //
        // SUBWORKFLOW: GENERATE SUPPLEMENTARY FILES FOR PRETEXT INGESTION
        //
        PRETEXT_ACCESSORY_FILES (
            ch_reference_tuple,
            ch_reference_sizes,
            ch_longread_reads,
            ch_teloseq,
            map_params["run_telomere"],
            map_params["run_gap"],
            map_params["run_coverage"],
            map_params["run_repeat_den"],
            map_params["run_busco"],
            map_params["run_pebble"],
            map_params["run_track_indexes"],
        )

        gaps_file       = PRETEXT_ACCESSORY_FILES.out.gap_file.map{ _meta, file -> file }.ifEmpty{ [] }
        cove_file       = PRETEXT_ACCESSORY_FILES.out.coverage_file.map{ _meta, file -> file }.ifEmpty{ [] }
        rept_file       = PRETEXT_ACCESSORY_FILES.out.repeat_file.map{ _meta, file -> file }.ifEmpty{ [] }
        telo_file       = PRETEXT_ACCESSORY_FILES.out.telo_windows.map{ _meta, file -> file }.collect().ifEmpty{ [] }
    }


    //
    // SUBWORKFLOW: MAP CRAM IF READS NOT ALREADY MAPPED
    //
    CRAM_MAP_ILLUMINA_HIC (
        ch_reference_tuple,
        ch_hic_reads,
        val_selected_aligner,
        val_cram_chunk_size
    )

    mapped_bam          = pre_mapped_bam.mix( CRAM_MAP_ILLUMINA_HIC.out.bam )


    //
    // MODULE: IF A PRE-MAPPED BAM HAS BEEN PROVIDED THEN INDEX IT
    //
    SAMTOOLS_INDEX (
        pre_mapped_bam
    )
    mapped_bam_bai      = SAMTOOLS_INDEX.out.index.mix( CRAM_MAP_ILLUMINA_HIC.out.bam_index )


    //
    // MODULE: SAMTOOLS FLAGSTAT TO COLLECT SOME MAPPING STATISTICS
    //
    mapped_bam_input      = mapped_bam.combine(mapped_bam_bai, by: 0)
    SAMTOOLS_FLAGSTAT (
        mapped_bam_input.filter { _meta, _bam, _csi -> val_mapping_statistics }
    )

    //
    // LOGIC: IF params.snapshot_order IS PROVIDED, USE IT TO ORDER SNAPSHOTS
    //        OTHERWISE, THE MODULE SHOULD STILL RUN WITHOUT ORDERING AND
    //        PRODUCE AN EMPTY CHANNEL. ONLY NEEDED FOR STDRD
    //
    ch_snapshot_custom_order = ch_reference_tuple
        .combine(ch_snapshot_order)
        .map { meta, _fasta, order_file -> [meta, order_file] }


    //
    // SUBWORKFLOW: MAP THE PRETEXT FILE AND TAKE SNAPSHOT
    //              STNDRD IS THE ONLY VARIANT WE ARE ANNOTATING WITH OTHER PARAMS
    //              DUE TO HOW RESOURCE INTENSIVE THEY SNAPSHOT IS WITH HIGHER RESOLUTION
    //              AND WE ONLY NEED 1 JUICER MAP
    //
    for_mapping = mapped_bam
        .combine(ch_reference_sizes, by: 0)

    // Turn off the tools which need pairs files
    CREATE_MAPS_STDRD (
        for_mapping.map { meta, bam, _sizes -> tuple(meta, bam) },
        for_mapping.map { meta, _bam, sizes -> tuple(meta, sizes) },
        ch_snapshot_custom_order,
        true,                       // Pretext generation, always true
        val_snapshot_generation,
        false,  //true,             // Cooler map generation, which we won't be using
        false,  // Juicer generation, optional need for genomenotes
        []                          // Cooler cload parameters
    )


    //
    // MODULE: ANNOTATE THE PRETEXT SNAPSHOT FILE WITH SCAFFOLD NAMES AND SIZES
    //
    filtered_sizes = ch_reference_sizes.filter { _meta, _file -> val_snapshot_annotate }

    PRETEXTANNOTATE(
        filtered_sizes,
        CREATE_MAPS_STDRD.out.pretext_png
    )


    //
    // SUBWORKFLOW: MAP THE PRETEXT FILE
    //
    CREATE_MAPS_HIRES (
        mapped_bam.filter{ map_params["run_hires"] },
        [[:],[]],
        channel.of([[:],[]]),
        true,
        false,
        false,
        false,
        []
    )

    //
    // SUBWORKFLOW: MAP THE PRETEXT FILE
    //              IF val_run_ultra IS "true" CALCULATE WHETHER THE REF IS > 4GB AND MAP ULTRA
    //              IF val_run_ultra IS "force" MAP ULTRA
    //
    def ultra_input     = mapped_bam
        .combine(ch_reference_tuple)
        .filter { _mapped_meta, _bam, _ref_meta, ref_fasta ->
            map_params["run_ultra"] == "force" || (map_params["run_ultra"] == "yes" && ref_fasta.size() > 4.GB)
        }
        .map { mapped_meta, bam, _ref_meta, _ref_fasta ->
            [mapped_meta, bam]
        }

    CREATE_MAPS_ULTRA (
        ultra_input,
        [[:],[]],
        channel.of([[:],[]]),
        true,
        false,
        false,
        false,
        []
    )


    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //          - ADAPTED FROM TREEVAL
    //
    PRETEXT_INGEST_SNDRD (
        CREATE_MAPS_STDRD.out.pretext.filter { map_params["run_generate_tracks"] },
        gaps_file,
        cove_file,
        telo_file,
        rept_file,
        map_params["run_split_telomere"]
    )


    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //          - ADAPTED FROM TREEVAL
    //
    PRETEXT_INGEST_HIRES (
        CREATE_MAPS_HIRES.out.pretext.filter { map_params["run_generate_tracks"] },
        gaps_file,
        cove_file,
        telo_file,
        rept_file,
        map_params["run_split_telomere"]
    )


    //
    // MODULE: INGEST ACCESSORY FILES INTO PRETEXT BY DEFAULT
    //          - ADAPTED FROM TREEVAL
    //
    PRETEXT_INGEST_ULTRA (
        CREATE_MAPS_ULTRA.out.pretext.filter { map_params["run_generate_tracks"] },
        gaps_file,
        cove_file,
        telo_file,
        rept_file,
        map_params["run_split_telomere"]
    )

}
