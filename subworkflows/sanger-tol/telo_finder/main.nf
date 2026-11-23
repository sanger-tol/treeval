//
// MODULE IMPORT BLOCK
//
include { TELOMERE_FINDTELOMERE } from '../../../modules/sanger-tol/telomere/findtelomere/main'
include { HTSLIB_BGZIPTABIX     } from '../../../modules/nf-core/htslib/bgziptabix/main'


workflow TELO_FINDER {

    take:
    ch_reference           // channel: [ val(meta), path(fasta) ]
    ch_telomereseq         // channel: [ val(meta), val(telomere_motif) ]
    val_zip_bed            // boolean — bgzip+tabix BED/windows (not the *.telomere TSV)

    main:

    ch_joined = ch_reference
        .combine(ch_telomereseq, by: 0)
        .map { meta, reference, telomereseq -> tuple(meta, reference, telomereseq) }

    TELOMERE_FINDTELOMERE(ch_joined)

    /*
     * Zip set mirrors FindTelomereWindows outputs from sanger-tol/telomere:
     *   always:  *.telomere.bed, *.windows
     *   --split: *.fwd/.rev.telomere.bed, *.fwd/.rev.windows
     * Optional channels emit nothing when split is false.
     * Never zip the *.telomere TSV (not BED / not windows).
     */
    ch_for_zip_raw = TELOMERE_FINDTELOMERE.out.telomere_bed
        .mix(TELOMERE_FINDTELOMERE.out.windows_all)
        .mix(TELOMERE_FINDTELOMERE.out.telomere_bed_fwd)
        .mix(TELOMERE_FINDTELOMERE.out.telomere_bed_rev)
        .mix(TELOMERE_FINDTELOMERE.out.windows_fwd)
        .mix(TELOMERE_FINDTELOMERE.out.windows_rev)

    if (val_zip_bed) {
        /*
         * HTSLIB_BGZIPTABIX expects
         *   [ val(meta), path(infile), path(infile_tbi), path(regions) ]
         * per element. Optional globs may yield a list of paths — expand those.
         */
        ch_for_zip = ch_for_zip_raw
            .flatMap { meta, item ->
                def items = (item instanceof List) ? item : [ item ]
                items.collect { file -> tuple(meta, file, [], []) }
            }

        HTSLIB_BGZIPTABIX(
            ch_for_zip,
            'compress',
            true,
            ''
        )

        ch_gz_index = HTSLIB_BGZIPTABIX.out.output
            .combine(HTSLIB_BGZIPTABIX.out.index)
            .filter { _meta, gz, _meta2, idx -> idx.name.startsWith(gz.name) }
            .map { meta, gz, _meta2, idx -> tuple(meta, gz, idx) }

    } else {
        ch_gz_index = channel.empty()
    }

    emit:
    telomere         = TELOMERE_FINDTELOMERE.out.telomere         // channel: [ val(meta), path(*.telomere) ]
    telomere_bed     = TELOMERE_FINDTELOMERE.out.telomere_bed     // channel: [ val(meta), path(*.telomere.bed) ]
    telomere_bed_fwd = TELOMERE_FINDTELOMERE.out.telomere_bed_fwd // channel: [ val(meta), path(*.fwd.telomere.bed) ]
    telomere_bed_rev = TELOMERE_FINDTELOMERE.out.telomere_bed_rev // channel: [ val(meta), path(*.rev.telomere.bed) ]
    windows_all      = TELOMERE_FINDTELOMERE.out.windows_all      // channel: [ val(meta), path(*.windows) ]
    windows_fwd      = TELOMERE_FINDTELOMERE.out.windows_fwd      // channel: [ val(meta), path(*.fwd.windows) ]
    windows_rev      = TELOMERE_FINDTELOMERE.out.windows_rev      // channel: [ val(meta), path(*.rev.windows) ]
    gz_index         = ch_gz_index                                 // channel: [ val(meta), path(*.gz), path(*.{tbi,csi}) ]

}
