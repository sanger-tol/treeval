class TreeValProject {

    //
    // Generate a summary file of input file data (size of files)
    // Still in testing
    //

    public static void summary(workflow, rf_data, pb_data, cm_data, summary_params, projectDir, params) {

        def input_data = [:]
        input_data['version']           = NfcoreTemplate.version( workflow )
        input_data['runName']           = workflow.runName
        input_data['session_id']        = workflow.sessionId
        input_data['duration']          = workflow.duration
        input_data['Date Started']      = workflow.start
        input_data['Date Completed']    = workflow.complete

        input_data['sample_id']         = rf_data.map{ it[0].id }
        input_data['taxonomic_class']   = rf_data.map{ it[0].ln }
        input_data['ticket_type']       = rf_data.map{ it[0].tk }
        input_data['input_asm']         = rf_data
        input_data['input_pacbio']      = ( pb_data ?: 'None'   )
        input_data['input_cram']        = ( cm_data ?: 'None'   )

        if (workflow.success) {
            def time = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')

            def output_directory = new File("${params.outdir}/pipeline_info/")
            if (!output_directory.exists()) {
                output_directory.mkdirs()
            }

            def output_hf = new File(output_directory, "input_data_${time}.txt")
            output_hf.withWriter { w -> w << input_data }

            new File( output_directory, "TreeVal_run_context_${time}.txt" ).withWriter { w ->
                ["---INPUT_DATA---\n",
                "${params.outdir}/pipeline_info/pipeline_execution_${time}.txt",
                "---RESOURCE_STATS---\n",
                "${params.outdir}/pipeline_info/input_data_${time}.txt"]
                    .each { f ->
                        new File( f ).withReader { r ->
                            w << r << '\n'
                        }
                    }
            }
        }
    }

    // Should generate a file looking like:
    // // // // // // //
    // ---INPUT_DATA---
    // Version:     {workflow.version}
    // runName:     {workflow.unName}
    // duration:    {workflow.duration}
    // input_asm:   [ [id, sz] file]
    // pacbio_bam:  [ [id, sz] file]
    // cram_bam:    [ [id, sz] file]
    // ---RESOURCE_STATS---
    // process  cpu cpu_usage   mem mem_usage   peak_usage
    // MINIMAP2_ALIGN   16  1590%   60GB    50% 30GB
    // ....

}
