include { PRETEXT_GRAPH                                 } from '../../modules/local/pretext/graph/main'

workflow PRETEXT_INGESTION {
    take:
    pretext_file        // Channel: tuple([sample_id], file)
    gap_file            // Channel: tuple([sample_id], file)
    coverage_file       // Channel: tuple([sample_id], file)
    cov_log_file        // Channel: tuple([sample_id], file)
    telomere_file       // Channel: tuple([sample_id], file)
    repeat_cov_file     // Channel: tuple([sample_id], file)


    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: GAP OR TELOMERE FILES CAN SOMETIMES BE EMPTY
    //          CHECK IF EMPTY AND ASSIGN APPROPRIATE BRANCHING
    //
    gap_file
        .map { meta, gap_file ->
            tuple( [    id: meta.id,
                        sz: gap_file.size().toInteger(),
                        ft: 'gap' ],
                        gap_file
            )
        }
        .set { ch_gap }

    telomere_file
        .map { meta, telo_file ->
            tuple( [    id: meta.id,
                        sz: telo_file.size().toInteger(),
                        ft: 'telomere' ],
                        telo_file
            )
        }
        .set { ch_telomere }

    //
    // MODULE: PRETEXT GRAPH INGESTS THE OTHER TWO FILES DIRECTLY INTO THE PRETEXT
    //          RUNNING AS IT'S OWN SUB IN ORDER TO NOT SLOW DOWN HIC_MAPPING ANY FURTHER
    //

    PRETEXT_GRAPH (
        pretext_file,
        ch_gap,
        coverage_file,
        cov_log_file,
        ch_telomere,
        repeat_cov_file
    )
    ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )

    emit:
    pretext_file        = PRETEXT_GRAPH.out.pretext
    versions            = ch_versions.ifEmpty(null)
}
