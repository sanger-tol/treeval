#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/treeval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/treeval
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_treeval_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_treeval_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { TREEVAL           } from './workflows/treeval'

//
// WORKFLOW: RUN MAIN PIPELINE GENERATING ALL OUTPUT
//
// workflow SANGERTOL_TREEVAL {
//     take:
//     assembly_id,
//     reference,
//     map_order,
//     assem_reads,
//     kmer_prof_file,
//     hic_reads,
//     supp_reads,
//     align_genesets,
//     synteny_paths,
//     intron_size,
//     teloseq,
//     lineageinfo,
//     lineagespath

//     main:

//     TREEVAL (
//         assembly_id,
//         reference,
//         map_order,
//         assem_reads,
//         kmer_prof_file,
//         hic_reads,
//         supp_reads,
//         align_genesets,
//         synteny_paths,
//         intron_size,
//         teloseq,
//         lineageinfo,
//         lineagespath
//     )
// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    params.mode    = params.mode ?: "FULL"

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.mode
    )

    //
    // WORKFLOW: Run main workflow
    //
    TREEVAL (
        PIPELINE_INITIALISATION.out.assembly_id,
        PIPELINE_INITIALISATION.out.reference,
        PIPELINE_INITIALISATION.out.map_order,
        PIPELINE_INITIALISATION.out.assem_reads,
        PIPELINE_INITIALISATION.out.kmer_prof_file,
        PIPELINE_INITIALISATION.out.hic_reads,
        PIPELINE_INITIALISATION.out.supp_reads,
        PIPELINE_INITIALISATION.out.align_genesets,
        PIPELINE_INITIALISATION.out.synteny_paths,
        PIPELINE_INITIALISATION.out.intron_size,
        PIPELINE_INITIALISATION.out.teloseq,
        PIPELINE_INITIALISATION.out.lineageinfo,
        PIPELINE_INITIALISATION.out.lineagespath
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
