/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// IMPORT: SUBWORKFLOWS CALLED BY THE MAIN
//
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

workflow TREEVAL {
    take:
    assembly_id     // channel:
    reference       // channel:
    map_order       // channel: hic mapping order (from yaml)
    assem_reads     // channel: path to longreads directory (from yaml)
    kmer_prof_file  // channel:
    hic_reads       // channel: path to hic reads directory (from yaml)
    supp_reads      // channel:
    align_genesets  // channel: paths to genesets in from yaml
    synteny_paths   // channel: path to syntenic genomes (from yaml)
    intron_size     // channel:
    teloseq         // channel: telomere motif sequence (from yaml)
    lineageinfo     // channel:
    lineagespath    // channel:

    main:
    //
    // PRE-PIPELINE CHANNEL SETTING - channel setting for required files
    //
    ch_versions         = Channel.empty()

    params.steps    = params.steps ?: 'NONE'
    exclude_steps_list = params.steps.length() > 1 ? params.steps.split(',').collect { it.trim() } : params.steps

    all_steps_list       = ["insilico_digest", "gene_alignment", "repeat_density", "gap_finder", "selfcomp", "synteny", "read_coverage", "telo_finder", "busco", "kmer", "hic_mapping", "NONE"]

    jbrowse_exclude_list = ["insilico_digest", "gene_alignment", "selfcomp", "synteny", "busco", "kmer"]
    rapid_exclude_list = ["gene_alignment", "repeat_density", "gap_finder", "read_coverage", "telo_finder", "hic_mapping"]
    rapid_tol_exclude_list = ["gene_alignment", "repeat_density", "gap_finder", "read_coverage", "telo_finder", "kmer",  "hic_mapping"]

    //
    // Add exclude determined by run mode (JBROWSE, RAPID, RAPID_TOL)
    //
    if (params.mode == "JBROWSE") {
        exclude_workflow_steps = (jbrowse_exclude_list + exclude_steps_list).unique()
    } else if (params.mode == "RAPID") {
        exclude_workflow_steps = (rapid_exclude_list + exclude_steps_list).unique()
    } else if (params.mode == "RAPID_TOL") {
        exclude_workflow_steps = (rapid_tol_exclude_list + exclude_steps_list).unique()
    } else {
        exclude_workflow_steps = exclude_steps_list
    }

    if (!all_steps_list.containsAll(exclude_workflow_steps)) {
        log.error "There is an extra argument given on Command Line (--steps): ${exclude_workflow_steps - all_steps_list}"
        log.error "Valid options are: ${all_steps_list.join(", ")}"
    }

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
    // SUBWORKFLOW: Takes input fasta file and sample ID to generate a my.genome file
    //
    GENERATE_GENOME (
        reference,
        map_order
    )
    ch_versions     = ch_versions.mix( GENERATE_GENOME.out.versions )


    //
    // SUBWORKFLOW: Takes reference, channel of enzymes, my.genome, assembly_id and as file to generate
    //              file with enzymatic digest sites.
    //
    if ( !(exclude_workflow_steps?.contains("insilico_digest"))) {
        ch_enzyme       = Channel.of( "bspq1","bsss1","DLE1" )

        INSILICO_DIGEST (
            GENERATE_GENOME.out.dot_genome,
            reference,
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
    if ( !(exclude_workflow_steps?.contains("gene_alignment"))) {
        GENE_ALIGNMENT (
            GENERATE_GENOME.out.dot_genome,
            reference,
            GENERATE_GENOME.out.ref_index,
            align_genesets,
            intron_size,
            gene_alignment_asfiles
        )
        ch_versions     = ch_versions.mix(GENE_ALIGNMENT.out.versions)
    }


    //
    // SUBWORKFLOW: GENERATES A BIGWIG FOR A REPEAT DENSITY TRACK
    //
    if ( !(exclude_workflow_steps?.contains("repeat_density"))) {
        REPEAT_DENSITY (
            reference,
            GENERATE_GENOME.out.dot_genome
        )
        ch_versions         = ch_versions.mix( REPEAT_DENSITY.out.versions )

        ch_repeat_density   = REPEAT_DENSITY.out.repeat_density
    } else {
        ch_repeat_density   = [[],[]]
    }


    //
    // SUBWORKFLOW: GENERATES A GAP.BED FILE TO ID THE LOCATIONS OF GAPS
    //
    if ( !(exclude_workflow_steps?.contains("gap_finder"))) {
        GAP_FINDER (
            reference
        )
        ch_versions         = ch_versions.mix( GAP_FINDER.out.versions )

        ch_gap_file         = GAP_FINDER.out.gap_file
    } else {
        ch_gap_file         = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: Takes reference file, .genome file, mummer variables, motif length variable and as
    //              file to generate a file containing sites of self-complementary sequnce.
    //
    if ( !(exclude_workflow_steps?.contains("selfcomp"))) {
        SELFCOMP (
            reference,
            GENERATE_GENOME.out.dot_genome,
            selfcomp_asfile
        )
        ch_versions         = ch_versions.mix( SELFCOMP.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reference, the directory of syntenic genomes and order/clade of sequence
    //              and generated a file of syntenic blocks.
    //
    if ( !(exclude_workflow_steps?.contains("synteny"))) {
        SYNTENY (
            reference,
            synteny_paths
        )
        ch_versions         = ch_versions.mix( SYNTENY.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reference, pacbio reads
    //
    if ( !(exclude_workflow_steps?.contains("read_coverage"))) {
        READ_COVERAGE (
            reference,
            GENERATE_GENOME.out.dot_genome,
            assem_reads
        )
        ch_versions         = ch_versions.mix( READ_COVERAGE.out.versions )
        ch_coverage_bg_norm = READ_COVERAGE.out.ch_covbw_nor
        ch_coverage_bg_avg  = READ_COVERAGE.out.ch_covbw_avg
    } else {
        ch_coverage_bg_avg  = Channel.of([[],[]])
        ch_coverage_bg_norm = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: GENERATE TELOMERE WINDOW FILES WITH PACBIO READS AND REFERENCE
    //
    if ( !(exclude_workflow_steps?.contains("telo_finder"))) {
        TELO_FINDER (   reference,
                        teloseq
        )
        ch_versions         = ch_versions.mix( TELO_FINDER.out.versions )
        ch_telo_bedgraph    = TELO_FINDER.out.bedgraph_file
    } else {
        ch_telo_bedgraph    = Channel.of([[],[]])
    }


    //
    // SUBWORKFLOW: GENERATE BUSCO ANNOTATION FOR ANCESTRAL UNITS
    //
    if ( !(exclude_workflow_steps?.contains("busco"))) {
        BUSCO_ANNOTATION (
            GENERATE_GENOME.out.dot_genome,
            reference,
            lineageinfo,
            lineagespath,
            buscogene_asfile,
            ancestral_table
        )
        ch_versions         = ch_versions.mix( BUSCO_ANNOTATION.out.versions )
    }


    //
    // SUBWORKFLOW: Takes reads and assembly, produces kmer plot
    //
    if ( !(exclude_workflow_steps?.contains("kmer"))) {
        KMER (
            reference,
            assem_reads
        )
        ch_versions         = ch_versions.mix( KMER.out.versions )
    }


    //
    // SUBWORKFLOW: GENERATE HIC MAPPING TO GENERATE PRETEXT FILES AND JUICEBOX
    //
    if ( !(exclude_workflow_steps?.contains("hic_mapping"))) {
        HIC_MAPPING (
            reference,
            GENERATE_GENOME.out.ref_index,
            GENERATE_GENOME.out.dot_genome,
            hic_reads,
            assembly_id,
            ch_gap_file,
            ch_coverage_bg_norm,
            ch_coverage_bg_avg,
            ch_telo_bedgraph,
            ch_repeat_density,
            params.mode
        )
        ch_versions         = ch_versions.mix( HIC_MAPPING.out.versions )
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
