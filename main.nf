#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/treeval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/treeval
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise( workflow, params, log )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TREEVAL           } from './workflows/treeval'
include { TREEVAL_RAPID     } from './workflows/treeval_rapid'
include { TREEVAL_RAPID_TOL } from './workflows/treeval_rapid_tol'

//
// WORKFLOW: RUN MAIN PIPELINE GENERATING ALL OUTPUT
//
workflow SANGERTOL_TREEVAL {
        TREEVAL ()
}

//
// WORKFLOW: RUN TRUNCATED PIPELINE TO PRODUCE CONTACT MAPS AND PRETEXT ACCESSORIES
//
workflow SANGERTOL_TREEVAL_RAPID {
        TREEVAL_RAPID ()
}

//
// WORKFLOW: RUN TRUNCATED PIPELINE, CONTAINS WORKFLOWS INTERNAL TO SANGERTOL
//
workflow SANGERTOL_TREEVAL_RAPID_TOL {
        TREEVAL_RAPID_TOL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOWS: Execute named workflow for the pipeline
//
workflow {
        SANGERTOL_TREEVAL ()
}

workflow RAPID {
        SANGERTOL_TREEVAL_RAPID ()
}

workflow RAPID_TOL {
        SANGERTOL_TREEVAL_RAPID_TOL ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
