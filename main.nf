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
include { TREEVAL_JBROWSE   } from './workflows/treeval_jbrowse'

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

//
// WORKFLOW: RUN ONLY THE SUBWORKFLOWS REQUIRED FOR JBROWSE UPLOAD
//              - THIS IS TO COMPLEMENT A NEW PROCESS WHERE MAJORITY OF TICKETS WILL BE RC
//                  AND GET REQUESTED FOR FULL
//
workflow SANGERTOL_TREEVAL_JBROWSE {
    TREEVAL_JBROWSE ()
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

workflow JBROWSE {
    SANGERTOL_TREEVAL_JBROWSE ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
