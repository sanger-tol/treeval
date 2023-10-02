include { PRETEXT_GRAPH } from '../../modules/local/pretext_graph'

workflow PRETEXT_INGESTION {
    take:
    pretext_file        // tuple([sample_id], file)
    gap_file            // tuple([sample_id], file)
    telomere_file       // tuple([sample_id], file)
    coverage_file       // tuple([sample_id], file)
    repeat_cov_file     // tuple([sample_id], file)

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: PRETEXT GRAPH INGESTS THE OTHER FOUR FILES DIRECTLY INTO THE PRETEXT
    //          RUNNING AS IT'S OWN SUB IN ORDER TO NOT SLOW DOWN HIC_MAPPING ANY FURTHER
    //Pre

    PRETEXT_GRAPH (
        pretext_file,
        gap_file,
        coverage_file,
        telomere_file,
        repeat_cov_file
    )
    ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )


    emit:
    pretext =   PRETEXT_GRAPH.out.pretext

}
