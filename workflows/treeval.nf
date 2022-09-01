/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTreeval.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GENERATE_GENOME   } from '../subworkflows/local/generate_genome'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow TREEVAL {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: reads the yaml and pushing out into a channel per yaml field
    //
    INPUT_READ ( params.input )

    //
    // SUBWORKFLOW: Takes input fasta file and sample ID to generate a my.genome file
    //    
    GENERATE_GENOME ( INPUT_READ.out.assembly_id, INPUT_READ.out.reference )
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    // USE GENERATE_GENOME.out.REFERENCE_TUPLE  // channel [[meta.id = sample], file(reference file)]
    // USE GENERATE_GENOME.out.dot_genome       // channel [[meta.id = sample], file(*.genome)]

    //
    //SUBWORKFLOW: 
    //
    //INSILICO_DIGEST ( INPUT_READ.out.sample_id,
    //                  GENERATE_GENOME.out.dot_genome,
    //                  GENERATE_GENOME.out.reference_tuple )
    //ch_versions = ch_versions.mix(INSILICO_DIGEST.out.versions)

    //
    //SUBWORKFLOW: Takes input fasta to generate BB files containing alignment data
    //
    //GENE_ALIGNMENT ( GENERATE_GENOME.out.dot_genome,
    //                 GENERATE_GENOME.out.reference_tuple,
    //                 INPUT_READ.out.assembly_classT,
    //                 INPUT_READ.out.align_data_dir,
    //                 INPUT_READ.out.align_geneset )
    //ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    //SUBWORKFLOW: 
    //
    //SELFCOMP ( GENERATE_GENOME.out.reference_tuple,
    //           GENERATE_GENOME.out.dot_genome,
    //           INPUT_READ.out.mummer_chunk,
    //           INPUT_READ.out.motif_len )
    //ch_versions = ch_versions.mix(SELFCOMP.out.versions)

    //
    //SUBWORKFLOW: 
    //
    //SYNTENY ( GENERATE_GENOME.out.reference_tuple )
    //ch_versions = ch_versions.mix(SYNTENY.out.versions)


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/