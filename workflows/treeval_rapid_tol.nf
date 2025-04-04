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
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_treeval_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TREEVAL_RAPID_TOL {

    main:
    ch_versions     = Channel.empty()

    exclude_workflow_steps  = params.exclude ? params.exclude.split(",") : "NONE"

    full_list       = ["insilico_digest", "gene_alignment", "repeat_density", "gap_finder", "selfcomp", "synteny", "read_coverage", "telo_finder", "busco", "kmer", "hic_mapping", "NONE"]

    if (!full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of --exclude: $exclude_workflow_steps\nMaster list is: $full_list"
    }

    params.entry    = 'RAPID_TOL'
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
    if ( !exclude_workflow_steps.contains("repeat_density")) {
        REPEAT_DENSITY (
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.dot_genome
        )
        ch_versions     = ch_versions.mix( REPEAT_DENSITY.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    if ( !exclude_workflow_steps.contains("gap_finder")) {
        GAP_FINDER (
            YAML_INPUT.out.reference_ch
        )
        ch_versions     = ch_versions.mix( GAP_FINDER.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH PACBIO READS AND REFERENCE
    //
    if ( !exclude_workflow_steps.contains("telo_finder")) {
        TELO_FINDER (   YAML_INPUT.out.reference_ch,
                        YAML_INPUT.out.teloseq
        )
        ch_versions     = ch_versions.mix( TELO_FINDER.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reference, pacbio reads
    //
    if ( !exclude_workflow_steps.contains("read_coverage")) {
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
    // SUBWORKFLOW: Takes reads and assembly, produces kmer plot
    //
    if ( !exclude_workflow_steps.contains("kmer")) {
        KMER (
            YAML_INPUT.out.reference_ch,
            YAML_INPUT.out.read_ch
        )
        ch_versions     = ch_versions.mix( KMER.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATE HIC MAPPING TO GENERATE PRETEXT FILES AND JUICEBOX
    //
    if ( !exclude_workflow_steps.contains("hic_mapping")) {
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
        hic_report      = HIC_MAPPING.out.ch_reporting
        ch_versions     = ch_versions.mix( HIC_MAPPING.out.versions )
    } else {
        hic_report = []
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'treeval_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
