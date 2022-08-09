//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

include { MAKECMAP_RENAMECMAPIDS } from '../modules/sanger-tol/nf-core-modules/makecmap/fa2cmapmulticolor/main'
include { MAKECMAP_RENAMECMAPIDS } from '../modules/sanger-tol/nf-core-modules/makecmap/renamecmapids/main'


workflow insilico_digest {

    main:
    ch_versions = Channel.empty()

    // Uncompress genome fasta file if required
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    
    
    emit:
    fasta    = REMOVE_MASKING.out.fasta  // path: genome.unmasked.fasta
    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

    enzyme = val(params.enzyme)

    //make cmap using fa2cmapmulticolor

    fa2cmapmulticolor(ch_fasta,enzyme)
    fa2cmapmulticolor.out.view()
}