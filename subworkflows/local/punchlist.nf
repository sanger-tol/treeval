include { PAFTOOLS_SAM2PAF      } from '../../modules/local/paftools_sam2paf'
include { PAF2BED               } from '../../modules/local/paf_to_bed12'

workflow PUNCHLIST {
    take:
    merged_bam

    main:
    ch_versions         = Channel.empty()

    //
    // MODULE: CONVERTS BAM INTO PAF FOR THE PUNCHLIST GENERATION
    //
    PAFTOOLS_SAM2PAF ( merged_bam )
    ch_versions     = ch_versions.mix(PAFTOOLS_SAM2PAF.out.versions)

    //
    // MODULE: GENERATES PUNCHLIST FROM PAF FILE
    //   Pass in prefix
    PAF2BED ( PAFTOOLS_SAM2PAF.out.paf )
    ch_versions     = ch_versions.mix(PAF2BED.out.versions)

    emit:
    punchlist   =  PAF2BED.out.bed
}
