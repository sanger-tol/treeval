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
include { READ_COVERAGE                                 } from '../subworkflows/local/read_coverage'
include { REPEAT_DENSITY                                } from '../subworkflows/local/repeat_density'
include { GAP_FINDER                                    } from '../subworkflows/local/gap_finder'
include { TELO_FINDER                                   } from '../subworkflows/local/telo_finder'
include { BUSCO_ANNOTATION                              } from '../subworkflows/local/busco_annotation'
include { HIC_MAPPING                                   } from '../subworkflows/local/hic_mapping'
include { PRETEXT_INGESTION as PRETEXT_INGEST_STANDRD   } from '../subworkflows/local/pretext_ingestion'
include { PRETEXT_INGESTION as PRETEXT_INGEST_HIGHRES   } from '../subworkflows/local/pretext_ingestion'
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

workflow TREEVAL {
    main:
    //
    // PRE-PIPELINE CHANNEL SETTING - channel setting for required files
    //
    ch_versions     = Channel.empty()

    exclude_workflow_steps  = params.steps ? params.steps.split(",") : "NONE"

    full_list       = ["insilico_digest", "gene_alignment", "repeat_density", "gap_finder", "selfcomp", "synteny", "read_coverage", "telo_finder", "busco", "kmer", "hic_mapping", "NONE"]

    if (!full_list.containsAll(exclude_workflow_steps)) {
        exit 1, "There is an extra argument given on Command Line: \n Check contents of --exclude: $exclude_workflow_steps\nMaster list is: $full_list"
    }

    params.entry    = 'FULL'
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
    if ( !exclude_workflow_steps.contains("insilico_digest")) {
        ch_enzyme       = Channel.of( "bspq1","bsss1","DLE1" )

        INSILICO_DIGEST (
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.reference_ch,
            ch_enzyme,
            digest_asfile
        )
        ch_versions     = ch_versions.mix( INSILICO_DIGEST.out.versions )
    }


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
    if ( !exclude_workflow_steps.contains("gene_alignment")) {
        GENE_ALIGNMENT (
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.ref_index,
            YAML_INPUT.out.align_genesets,
            YAML_INPUT.out.intron_size,
            gene_alignment_asfiles
        )
        ch_versions     = ch_versions.mix(GENE_ALIGNMENT.out.versions)
    }


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
    // SUBWORKFLOW: Takes reference file, .genome file, mummer variables, motif length variable and as
    //              file to generate a file containing sites of self-complementary sequnce.
    //
    if ( !exclude_workflow_steps.contains("selfcomp")) {
        SELFCOMP (
            YAML_INPUT.out.reference_ch,
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.motif_len,
            selfcomp_asfile
        )
        ch_versions     = ch_versions.mix( SELFCOMP.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reference, the directory of syntenic genomes and order/clade of sequence
    //              and generated a file of syntenic blocks.
    //
    if ( !exclude_workflow_steps.contains("synteny")) {
        SYNTENY (
            YAML_INPUT.out.reference_ch,
            YAML_INPUT.out.synteny_paths
        )
        ch_versions     = ch_versions.mix( SYNTENY.out.versions )
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
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH PACBIO READS AND REFERENCE
    //
    if ( !exclude_workflow_steps.contains("telo_finder")) {
        TELO_FINDER (   YAML_INPUT.out.reference_ch,
                        YAML_INPUT.out.teloseq
        )
        ch_versions     = ch_versions.mix( TELO_FINDER.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATE BUSCO ANNOTATION FOR ANCESTRAL UNITS
    //
    if ( !exclude_workflow_steps.contains("busco")) {
        BUSCO_ANNOTATION (
            GENERATE_GENOME.out.dot_genome,
            YAML_INPUT.out.reference_ch,
            YAML_INPUT.out.lineageinfo,
            YAML_INPUT.out.lineagespath,
            buscogene_asfile,
            ancestral_table
        )
        ch_versions = ch_versions.mix( BUSCO_ANNOTATION.out.versions )
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

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }

    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
