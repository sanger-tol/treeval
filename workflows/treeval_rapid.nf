/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTreeval.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// IMPORT: SUBWORKFLOWS CALLED BY THE MAIN
//
include { YAML_INPUT                                    } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME                               } from '../subworkflows/local/generate_genome'
include { REPEAT_DENSITY                                } from '../subworkflows/local/repeat_density'
include { GAP_FINDER                                    } from '../subworkflows/local/gap_finder'
include { READ_COVERAGE                                 } from '../subworkflows/local/read_coverage'
include { TELO_FINDER                                   } from '../subworkflows/local/telo_finder'
include { HIC_MAPPING                                   } from '../subworkflows/local/hic_mapping'
include { KMER                                          } from '../subworkflows/local/kmer'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// IMPORT: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TREEVAL_RAPID {

    main:
    ch_versions     = Channel.empty()

    exclude_workflow_steps = params.steps.length() > 1 ? params.steps.split(',').collect { it.trim() } : params.steps

    full_list       = ["insilico_digest", "gene_alignment", "repeat_density", "gap_finder", "selfcomp", "synteny", "read_coverage", "telo_finder", "busco", "kmer", "hic_mapping", "NONE"]

    if (!full_list.containsAll(exclude_workflow_steps)) {
        log.error "There is an extra argument given on Command Line (--steps): ${exclude_workflow_steps - full_list}"
        log.error "Valid options are: ${full_list.join(", ")}"
    }


    params.entry    = 'RAPID'
    input_ch        = Channel.fromPath(params.input, checkIfExists: true)


    //
    // SUBWORKFLOW: reads the yaml and pushing out into a channel per yaml field
    //
    YAML_INPUT (
        input_ch,
        params.entry
    )


    //
    // SUBWORKFLOW: Takes input fasta file and sample ID to generate a my.genome file
    //
    GENERATE_GENOME (
        YAML_INPUT.out.reference_ch,
        YAML_INPUT.out.map_order_ch
    )
    ch_versions     = ch_versions.mix( GENERATE_GENOME.out.versions )


    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    if ( !(exclude_workflow_steps?.contains("repeat_density"))) {
        REPEAT_DENSITY (
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.dot_genome
        )
        ch_versions     = ch_versions.mix( REPEAT_DENSITY.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    if ( !(exclude_workflow_steps?.contains("gap_finder"))) {
        GAP_FINDER (
            YAML_INPUT.out.reference_ch
        )
        ch_versions     = ch_versions.mix( GAP_FINDER.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH PACBIO READS AND REFERENCE
    //
    if ( !(exclude_workflow_steps?.contains("telo_finder"))) {
        TELO_FINDER (   YAML_INPUT.out.reference_ch,
                        YAML_INPUT.out.teloseq
        )
        ch_versions     = ch_versions.mix( TELO_FINDER.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reference, pacbio reads
    //
    if ( !(exclude_workflow_steps?.contains("read_coverage"))) {
        READ_COVERAGE (
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.read_ch
        )
        coverage_report = READ_COVERAGE.out.ch_reporting
        ch_versions     = ch_versions.mix( READ_COVERAGE.out.versions )
    } else {
        coverage_report = []
    }


    //
    // SUBWORKFLOW: GENERATE HIC MAPPING TO GENERATE PRETEXT FILES AND JUICEBOX
    //
    if ( !(exclude_workflow_steps?.contains("hic_mapping"))) {
        HIC_MAPPING (
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.ref_index,
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.hic_reads_ch,
            YAML_INPUT.out.assembly_id,
            GAP_FINDER.out.gap_file,
            READ_COVERAGE.out.ch_covbw_nor,
            READ_COVERAGE.out.ch_covbw_avg,
            TELO_FINDER.out.bedgraph_file,
            REPEAT_DENSITY.out.repeat_density,
            params.entry
        )
        hic_report = HIC_MAPPING.out.ch_reporting
        ch_versions     = ch_versions.mix( HIC_MAPPING.out.versions )
    } else {
        hic_report = []
    }


    //
    // SUBWORKFLOW: Collates version data from prior subworflows
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // LOGIC: GENERATE SOME CHANNELS FOR REPORTING
    //
    YAML_INPUT.out.reference_ch
        .combine( coverage_report )
        .combine( hic_report )
        .combine( CUSTOM_DUMPSOFTWAREVERSIONS.out.versions )
        .map { meta, reference, read_meta, read_files, hic_meta, hic_files, custom_file -> [
            rf_data: tuple(
                [   id: meta.id,
                    sz: file(reference).size(),
                    ln: meta.class,
                    tk: meta.project_id  ],
                reference
            ),
            sample_id: meta.id,
            pb_data: tuple( read_meta, read_files ),
            cm_data: tuple( hic_meta, hic_files ),
            custom: custom_file,
            ]
        }
        .set { collected_metrics_ch }

    collected_metrics_ch.map { metrics ->
        TreeValProject.summary( workflow, params, metrics, log )
    }

    emit:
    software_ch     = CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
    versions_ch     = CUSTOM_DUMPSOFTWAREVERSIONS.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// PIPELINE ENTRYPOINT SUBWORKFLOWS WILL USE THE IMPLICIT ONCOMPLETE BLOCK

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
