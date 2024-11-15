/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowTreeval.initialise(params, log)

// Check input path parameters to see if they exist
// params.input is the treeval yaml
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
include { INSILICO_DIGEST                               } from '../subworkflows/local/insilico_digest'
include { GENE_ALIGNMENT                                } from '../subworkflows/local/gene_alignment'
include { SELFCOMP                                      } from '../subworkflows/local/selfcomp'
include { SYNTENY                                       } from '../subworkflows/local/synteny'
include { BUSCO_ANNOTATION                              } from '../subworkflows/local/busco_annotation'
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

workflow TREEVAL_JBROWSE {
    main:
    //
    // PRE-PIPELINE CHANNEL SETTING - channel setting for required files
    //
    ch_versions     = Channel.empty()

    params.entry    = 'JBROWSE'
    input_ch        = Channel.fromPath(params.input, checkIfExists: true)

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

    Channel
        .fromPath( "${projectDir}/assets/busco_gene/busco.as", checkIfExists: true )
        .set { buscogene_asfile }

    Channel
        .fromPath( "${projectDir}/assets/busco_gene/lep_ancestral.tsv", checkIfExists: true )
        .set { ancestral_table }

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
    // SUBWORKFLOW: Takes reference, channel of enzymes, my.genome, assembly_id and as file to generate
    //              file with enzymatic digest sites.
    //
    ch_enzyme       = Channel.of( "bspq1","bsss1","DLE1" )

    INSILICO_DIGEST (
        GENERATE_GENOME.out.dot_genome,
        YAML_INPUT.out.reference_ch,
        ch_enzyme,
        digest_asfile
    )
    ch_versions     = ch_versions.mix( INSILICO_DIGEST.out.versions )

    //
    // SUBWORKFLOW: FOR SPLITTING THE REF GENOME INTO SCAFFOLD CHUNKS AND RUNNING SOME SUBWORKFLOWS
    //              ON THOSE CHUNKS
    //              THIS WILL BE REQUIRED FOR LARGER GENOMES EST > 6GB
    //
    // REFERENCE_GENOME_SPLIT --> SELFCOMP
    //                        --> GENE_ALIGNMENT
    //              BOTH WOULD REQUIRE A POST SUBWORKFLOW MERGE STEP TO MERGE TOGETHER THE SCAFFOLD
    //              BASED ALIGNMENTS/SELFCOMPS INTO A GENOME REPRESENTATIVE ONE.
    //              FOR GENE ALIGNMENT WOULD THIS REQUIRE A .GENOME FILE AND INDEX PER SCAFFOLD?

    //
    // SUBWORKFLOW: Takes input fasta to generate BB files containing alignment data
    //
    GENE_ALIGNMENT (
        GENERATE_GENOME.out.dot_genome,
        YAML_INPUT.out.reference_ch,
        GENERATE_GENOME.out.ref_index,
        YAML_INPUT.out.align_data_dir,
        YAML_INPUT.out.align_geneset,
        YAML_INPUT.out.align_common,
        YAML_INPUT.out.intron_size,
        gene_alignment_asfiles
    )
    ch_versions     = ch_versions.mix(GENE_ALIGNMENT.out.versions)

    //
    // SUBWORKFLOW: Takes reference file, .genome file, mummer variables, motif length variable and as
    //              file to generate a file containing sites of self-complementary sequnce.
    //
    SELFCOMP (
        YAML_INPUT.out.reference_ch,
        GENERATE_GENOME.out.dot_genome,
        YAML_INPUT.out.mummer_chunk,
        YAML_INPUT.out.motif_len,
        selfcomp_asfile
    )
    ch_versions     = ch_versions.mix( SELFCOMP.out.versions )

    //
    // SUBWORKFLOW: Takes reference, the directory of syntenic genomes and order/clade of sequence
    //              and generated a file of syntenic blocks.
    //
    SYNTENY (
        YAML_INPUT.out.reference_ch,
        YAML_INPUT.out.synteny_path
    )
    ch_versions     = ch_versions.mix( SYNTENY.out.versions )

    //
    // SUBWORKFLOW: GENERATE BUSCO ANNOTATION FOR ANCESTRAL UNITS
    //
    BUSCO_ANNOTATION (
        GENERATE_GENOME.out.dot_genome,
        YAML_INPUT.out.reference_ch,
        YAML_INPUT.out.lineageinfo,
        YAML_INPUT.out.lineagespath,
        buscogene_asfile,
        ancestral_table
    )
    ch_versions = ch_versions.mix( BUSCO_ANNOTATION.out.versions )

    //
    // SUBWORKFLOW: Takes reads and assembly, produces kmer plot
    //
    KMER (
        YAML_INPUT.out.reference_ch,
        YAML_INPUT.out.read_ch
    )
    ch_versions     = ch_versions.mix( KMER.out.versions )

    //
    // SUBWORKFLOW: Collates version data from prior subworflows
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
