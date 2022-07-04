//
// Synteny
//

include { MINIMAP2_ALIGN } from '../../modules/local/minimap2/align'

workflow SYNTENY {
    take:
    reads // channel: [ val(meta), [ datafile ] ]
    index // channel: /path/to/mmi
    fasta // channel: /path/to/fasta

    main:
    ch_versions = Channel.empty()

    // Check if species group has references for comparison

    //minimap2 -t 8 -x asm10 reference.fa my.fa > ref_my.paf
    // -t - number of threads
    // -x - preset mode
    // asm10 - Long assembly to reference mapping, Typically, the alignment will not extend to regions with 10% or higher sequence divergence. Only use this preset if the average divergence is far below 10%.

    // Align Fastq to reference
    MINIMAP2_ALIGN ( reads, fasta, index )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    emit:
    paf = MINIMAP2_ALIGN.out.paf

    versions = ch_versions
}

/* //
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}
