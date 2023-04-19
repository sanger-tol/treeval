include { PAFTOOLS_SAM2PAF      } from '../../modules/nf-core/paftools/sam2paf/main'
include { PAF2BED               } from '../../modules/local/paf_to_bed12'
//include { MV_TO_SANGER   } from '../../modules/local/mv_to_sanger'

workflow PUNCHLIST {
    take:
    reference_tuple
    dbVersion
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

    //
    // MODULE: TEMPORARY MODULE TO COPY THE PUNCHLIST DATA TO SANGER DIR
    //
    //  MV_TO_SANGER (  reference_tuple.map { it[0] },
    //                  dbVersion,
    //                  PAF2BED.out.punchlist )

    emit:
    punchlist   = PAF2BED.out.punchlist
    versions    = ch_versions.ifEmpty(null)
}
