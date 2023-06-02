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
// SUBWORKFLOW TREEVAL_FULL: Consisting of a mix of local and nf-core/modules
//
include { YAML_INPUT        } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME   } from '../subworkflows/local/generate_genome'
include { INSILICO_DIGEST   } from '../subworkflows/local/insilico_digest'
include { GENE_ALIGNMENT    } from '../subworkflows/local/gene_alignment'
include { SELFCOMP          } from '../subworkflows/local/selfcomp'
include { SYNTENY           } from '../subworkflows/local/synteny'
include { REPEAT_DENSITY    } from '../subworkflows/local/repeat_density'
include { GAP_FINDER        } from '../subworkflows/local/gap_finder'
include { LONGREAD_COVERAGE } from '../subworkflows/local/longread_coverage'
include { TELO_FINDER       } from '../subworkflows/local/telo_finder'
include { HIC_MAPPING       } from '../subworkflows/local/hic_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TREEVAL_RAPID {
    
    main:
    ch_versions = Channel.empty()

    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    //
    // SUBWORKFLOW: reads the yaml and pushing out into a channel per yaml field
    //
    YAML_INPUT ( input_ch )

    //
    // SUBWORKFLOW: Takes input fasta file and sample ID to generate a my.genome file
    //    
    GENERATE_GENOME ( YAML_INPUT.out.assembly_id,
                      YAML_INPUT.out.reference
    )
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    REPEAT_DENSITY ( GENERATE_GENOME.out.reference_tuple,
                     GENERATE_GENOME.out.dot_genome
    )
    ch_versions = ch_versions.mix(REPEAT_DENSITY.out.versions)

    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    GAP_FINDER ( GENERATE_GENOME.out.reference_tuple,
                 GENERATE_GENOME.out.max_scaff_size
    )
    ch_versions = ch_versions.mix(GAP_FINDER.out.versions)

    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH PACBIO READS AND REFERENCE
    //
    TELO_FINDER (       GENERATE_GENOME.out.reference_tuple,
                        YAML_INPUT.out.teloseq
    )
    ch_versions = ch_versions.mix(TELO_FINDER.out.versions) 

    //
    // SUBWORKFLOW: GENERATE HIC MAPPING TO GENERATE PRETEXT FILES AND JUICEBOX
    //
    HIC_MAPPING ( GENERATE_GENOME.out.reference_tuple,
                  GENERATE_GENOME.out.ref_index,
                  GENERATE_GENOME.out.dot_genome,
                  YAML_INPUT.out.hic_reads)
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    //
    // SUBWORKFLOW: Takes reference, pacbio reads 
    //
    LONGREAD_COVERAGE ( GENERATE_GENOME.out.reference_tuple,
                        GENERATE_GENOME.out.dot_genome,
                        YAML_INPUT.out.pacbio_reads,
                        YAML_INPUT.out.assembly_sizeClass
    )
    ch_versions = ch_versions.mix(LONGREAD_COVERAGE.out.versions)
 
    //
    // SUBWORKFLOW: Collates version data from prior subworflows
    //

    emit:
    ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
