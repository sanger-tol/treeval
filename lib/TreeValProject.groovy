class TreeValProject {

    //
    // Generate a summary file of input file data (size of files)
    // Still in testing
    //

    public static void summary(workflow, params) {

        def input_data = [:]
        input_data['version']           = NfcoreTemplate.version( workflow )
        input_data['runName']           = workflow.runName
        input_data['session_id']        = workflow.sessionId
        input_data['duration']          = workflow.duration
        input_data['DateStarted']       = workflow.start
        input_data['DateCompleted']     = workflow.complete

        input_data['rf_data']           = params.rf_data.value
        input_data['pb_data']           = 'None'
        input_data['cm_data']           = 'None'

        if (workflow.success) {
            def time = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')

            def output_directory = new File("${params.tracedir}/")
            if (!output_directory.exists()) {
                output_directory.mkdirs()
            }

            def output_hf = new File(output_directory, "input_data_${params.trace_timestamp}.txt")
            output_hf.write """\
                            ---RUN_DATA---
                            Pipeline_version:  ${input_data.version}
                            Pipeline_runname:  ${input_data.runName}
                            Pipeline_session:  ${input_data.session_id}
                            Pipeline_duration: ${input_data.duration}
                            Pipeline_datastrt: ${input_data.DateStarted}
                            Pipeline_datecomp: ${input_data.DateCompleted}
                            ---INPUT_DATA---
                            InputAssemblyData: ${input_data.rf_data}
                            Input_PacBio>Bam:  ${input_data.pb_data}
                            Input_HiCCram>Bam: ${input_data.cm_data}
                            ---RESOURCES---
                            """.stripIndent()

            def full_file = new File( output_directory, "TreeVal_run_context_${time}.txt" )
            def file_locs = ["${params.tracedir}/input_data_${time}.txt",
                                "${params.tracedir}/pipeline_execution_${params.trace_timestamp}.txt"]
            file_locs.each{ full_file.append( new File( it ).getText() ) }

        }
    }

    // Should generate a file looking like:
    // // // // // // //
    // ---INPUT_DATA---
    // Version:     {workflow.version}
    // runName:     {workflow.unName}
    // duration:    {workflow.duration}
    // ---INPUT_DATA---
    // input_asm:   [ [id, sz] file]
    // pacbio_bam:  [ [id, sz] file]
    // cram_bam:    [ [id, sz] file]
    // ---RESOURCE_STATS---
    // process  cpu cpu_usage   mem mem_usage   peak_usage
    // MINIMAP2_ALIGN   16  1590%   60GB    50% 30GB
    // ....
    // // // // // // //
}
