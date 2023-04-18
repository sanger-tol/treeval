/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTreeval.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { YAML_INPUT        } from '../subworkflows/local/yaml_input'
include { GENERATE_GENOME   } from '../subworkflows/local/generate_genome'
include { INSILICO_DIGEST   } from '../subworkflows/local/insilico_digest'
include { GENE_ALIGNMENT    } from '../subworkflows/local/gene_alignment'
include { SELFCOMP          } from '../subworkflows/local/selfcomp'
include { SYNTENY           } from '../subworkflows/local/synteny'
include { REPEAT_DENSITY    } from '../subworkflows/local/repeat_density'
include { GAP_FINDER        } from '../subworkflows/local/gap_finder'
// include { LONGREAD_COVERAGE } from '../subworkflows/local/longread_coverage'

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

workflow TREEVAL {

    //
    // PRE-PIPELINE CHANNEL SETTING - channel setting for required files
    //
    ch_versions = Channel.empty()

    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    Channel
        .fromPath( "${projectDir}/assets/gene_alignment/assm_*.as", checkIfExists: true)
        .map { it -> 
            tuple ([ type    :   it.toString().split('/')[-1].split('_')[-1].split('.as')[0] ],
                    file(it)
                )}
        .set { gene_alignment_asfiles }
    
    Channel
        .fromPath( "${projectDir}/assets/digest/digest.as", checkIfExists: true )
        .set { digest_asfile }

    Channel
        .fromPath( "${projectDir}/assets/self_comp/selfcomp.as", checkIfExists: true )
        .set { selfcomp_asfile }

    //
    // SUBWORKFLOW: reads the yaml and pushing out into a channel per yaml field
    //
    YAML_INPUT ( input_ch )

    //
    // SUBWORKFLOW: Takes input fasta file and sample ID to generate a my.genome file
    //    
    GENERATE_GENOME ( YAML_INPUT.out.assembly_id, YAML_INPUT.out.reference )
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)


    //
    // SUBWORKFLOW: Takes reference, channel of enzymes, my.genome, assembly_id and as file to generate
    //              file with enzymatic digest sites.
    //
    ch_enzyme = Channel.of( "bspq1","bsss1","DLE1" )
    INSILICO_DIGEST ( YAML_INPUT.out.assembly_id,
                      GENERATE_GENOME.out.dot_genome,
                      GENERATE_GENOME.out.reference_tuple,
                      ch_enzyme,
                      digest_asfile )
    ch_versions = ch_versions.mix(INSILICO_DIGEST.out.versions)

    //
    // SUBWORKFLOW: FOR SPLITTING THE REF GENOME INTO SCAFFOLD CHUNKS AND RUNNING SOME SUBWORKFLOWS
    //              ON THOSE CHUNKS
    //
    // REFERENCE_GENOME_SPLIT --> SELFCOMP
    //                        --> GENE_ALIGNMENT
    //              BOTH WOULD REQUIRE A POST SUBWORKFLOW MERGE STEP TO MERGE TOGETHER THE SCAFFOLD
    //              BASED ALIGNMENTS/SELFCOMPS INTO A GENOME REPRESENTATIVE ONE.
    //              FOR GENE ALIGNMENT WOULD THIS REQUIRE A .GENOME FILE AND INDEX PER SCAFFOLD?
 
    //
    // SUBWORKFLOW: Takes input fasta to generate BB files containing alignment data
    //
    GENE_ALIGNMENT ( GENERATE_GENOME.out.dot_genome,
                     GENERATE_GENOME.out.reference_tuple,
                     GENERATE_GENOME.out.ref_index,
                     YAML_INPUT.out.assembly_classT,
                     YAML_INPUT.out.align_data_dir,
                     YAML_INPUT.out.align_geneset,
                     YAML_INPUT.out.align_common,
                     YAML_INPUT.out.intron_size,
                     gene_alignment_asfiles,
                     YAML_INPUT.out.dbVersion )
    
    ch_versions = ch_versions.mix(GENERATE_GENOME.out.versions)

    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    REPEAT_DENSITY ( GENERATE_GENOME.out.reference_tuple,
                     GENERATE_GENOME.out.dot_genome )

    ch_versions = ch_versions.mix(REPEAT_DENSITY.out.versions)

    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    GAP_FINDER ( GENERATE_GENOME.out.reference_tuple )

    ch_versions = ch_versions.mix(GAP_FINDER.out.versions)

    //
    // SUBWORKFLOW: Takes reference file, .genome file, mummer variables, motif length variable and as
    //              file to generate a file containing sites of self-complementary sequnce.
    //
    SELFCOMP ( GENERATE_GENOME.out.reference_tuple,
               GENERATE_GENOME.out.dot_genome,
               YAML_INPUT.out.mummer_chunk,
               YAML_INPUT.out.motif_len,
               selfcomp_asfile )
    ch_versions = ch_versions.mix(SELFCOMP.out.versions)
 
    //
    // SUBWORKFLOW: Takes reference, the directory of syntenic genomes and order/clade of sequence
    //              and generated a file of syntenic blocks.
    //
    SYNTENY ( GENERATE_GENOME.out.reference_tuple, 
              YAML_INPUT.out.synteny_path,  
              YAML_INPUT.out.assembly_classT)
    ch_versions = ch_versions.mix(SYNTENY.out.versions)

    //
    // SUBWORKFLOW: 
    //
    // LONGREAD_COVERAGE (  GENERATE_GENOME.out.reference_tuple,
    //                      PACBIO.READ.DIRECTORY,
    //                      YAML_INPUT.out.sizeClass )

    //
    // SUBWORKFLOW: Collates version data from prior subworflows
    //
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
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
