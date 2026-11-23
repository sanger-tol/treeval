//
// GENERATE BED FILE OF GAPS AND LENGTH IN REFERENCE
//

include { SEQTK_CUTN                } from '../../../modules/nf-core/seqtk/cutn/main'
include { GAWK as GAWK_GAP_LENGTH   } from '../../../modules/nf-core/gawk/main'
include { HTSLIB_BGZIPTABIX          } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow GAP_FINDER {
    take:
    ch_reference     // Channel [ val(meta), path(fasta) ]
    val_run_bgzip    // val(boolean)

    main:

    //
    // MODULE: GENERATES A GAP SUMMARY FILE
    //
    SEQTK_CUTN (
        ch_reference
    )

    ch_reformat_gaps = channel.of('''\
        BEGIN { OFS = "\\t" } {
            print $0, sqrt(($3-$2)*($3-$2))
        }'''.stripIndent())
        .collectFile(name: "reformat_gaps.awk", cache: true)
        .collect()


    //
    // MODULE: ADD THE LENGTH OF GAP TO BED FILE - INPUT FOR PRETEXT MODULE
    //
    GAWK_GAP_LENGTH (
        SEQTK_CUTN.out.bed.map { meta, files -> tuple(meta, files, 'gap.bedgraph') },
        ch_reformat_gaps,
        false
    )


    //
    // MODULE: BGZIP AND TABIX THE GAP FILE
    //
    HTSLIB_BGZIPTABIX (
        SEQTK_CUTN.out.bed
            .map { meta, file -> tuple(meta, file, [], []) }
            .filter{ _meta, _file, _empty1, _empty2 -> val_run_bgzip},
        "compress",
        true,
        "bed"
    )

    htslib_bed_index = HTSLIB_BGZIPTABIX.out.output
        .combine(HTSLIB_BGZIPTABIX.out.index, by: 0)

    emit:
    gap_file        = GAWK_GAP_LENGTH.out.output
    gap_and_index   = htslib_bed_index
}
