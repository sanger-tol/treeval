include { PRETEXT_GRAPH                                 } from '../../modules/local/pretext_graph'

workflow PRETEXT_INGESTION {
    take:
    pretext_file        // tuple([sample_id], file)
    gap_file            // tuple([sample_id], file)
    coverage_file       // tuple([sample_id], file)
    telomere_file       // tuple([sample_id], file)
    repeat_cov_file     // tuple([sample_id], file)
    pretext_type        // var ( "hr || nr" )

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: GAP OR TELOMERE FILES CAN SOMETIMES BE EMPTY
    //          CHECK IF EMPTY AND ASSIGN APPROPRIATE BRANCHING
    //
    gap_file
        .map { meta, gap_file ->
            tuple( [    id: meta.id,
                        sz: gap_file.size(),
                        ft: 'gap' ],
                        gap_file
            )
        }
        .set { ch_gap }

    telomere_file
        .map { meta, telomere_file ->
            tuple( [    id: meta.id,
                        sz: telomere_file.size(),
                        ft: 'telomere' ],
                        telomere_file
            )
        }
        .set { ch_telomere }

    //
    // MODULE: PRETEXT GRAPH INGESTS THE OTHER TWO FILES DIRECTLY INTO THE PRETEXT
    //          RUNNING AS IT'S OWN SUB IN ORDER TO NOT SLOW DOWN HIC_MAPPING ANY FURTHER
    //

    PRETEXT_GRAPH (
        pretext_file,
        coverage_file,
        repeat_cov_file,
        [[],[]],
        ch_gap,
        ch_telomere,
        pretext_type
    )
    ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )

}
