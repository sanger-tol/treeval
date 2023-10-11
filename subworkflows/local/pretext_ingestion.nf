include { PRETEXT_GRAPH                                 } from '../../modules/local/pretext_graph'
include { PRETEXT_GRAPH_BG as PRETEXT_GRAPH_BG_GAP      } from '../../modules/local/pretext_graph_bg'
include { PRETEXT_GRAPH_BG as PRETEXT_GRAPH_BG_TELOMERE } from '../../modules/local/pretext_graph_bg'
include { PRETEXT_GRAPH_BG as PRETEXT_GRAPH_BG_BOTH_1   } from '../../modules/local/pretext_graph_bg'
include { PRETEXT_GRAPH_BG as PRETEXT_GRAPH_BG_BOTH_2   } from '../../modules/local/pretext_graph_bg'


workflow PRETEXT_INGESTION {
    take:
    pretext_file        // tuple([sample_id], file)
    gap_file            // tuple([sample_id], file)
    coverage_file       // tuple([sample_id], file)
    telomere_file       // tuple([sample_id], file)
    repeat_cov_file     // tuple([sample_id], file)

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
        .branch {
            valid:      it[0].sz >= 1
            invalid:    it[0].sz < 1
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
        .branch {
            valid:      it[0].sz >= 1
            invalid:    it[0].sz < 1
        }
        .set { ch_telomere }

    //
    // MODULE: PRETEXT GRAPH INGESTS THE OTHER TWO FILES DIRECTLY INTO THE PRETEXT
    //          RUNNING AS IT'S OWN SUB IN ORDER TO NOT SLOW DOWN HIC_MAPPING ANY FURTHER
    //

    PRETEXT_GRAPH (
        pretext_file,
        coverage_file,
        repeat_cov_file
    )
    ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )

    //
    // MODULE: PRETEXT GRAPH GB INGESTS
    //

    if ( ch_gap.valid.ifEmpty("YES") != "YES" && ch_telomere.valid.ifEmpty("YES") == "YES" ) {
        PRETEXT_GRAPH_BG_GAP (
            PRETEXT_GRAPH.out.pretext,
            ch_gap.valid
        )
        ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )
    } else if ( ch_gap.valid.ifEmpty("YES") == "YES" && ch_telomere.valid.ifEmpty("YES") != "YES" ) {
        PRETEXT_GRAPH_BG_TELOMERE (
            PRETEXT_GRAPH.out.pretext,
            ch_telomere.valid
        )
        ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )
    } else if ( ch_gap.valid.ifEmpty("YES") != "YES" && ch_telomere.valid.ifEmpty("YES") != "YES" ) {
        PRETEXT_GRAPH_BG_BOTH_1 (
            PRETEXT_GRAPH.out.pretext,
            ch_gap.valid
        )
        ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )

        PRETEXT_GRAPH_BG_BOTH_2 (
            PRETEXT_GRAPH_BG_BOTH_1.out.pretext,
            ch_telomere.valid
        )
        ch_versions         = ch_versions.mix( PRETEXT_GRAPH.out.versions )
    } else {
        println("BOTH FILES PRETEXT ACCESSORY FILES ARE EMPTY - only coverage and repeat density will be ingested")
    }
}
