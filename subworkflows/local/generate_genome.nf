include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { GENERATE_GENOME       } from '../../modules/local/genome_generator'

workflow GENERATE_GENOME {
    ch_versions     = Channel.empty()

    SAMTOOLS_FAIDX ( [[params.assembly.sample], params.reference] )
    ch_versions       = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GENERATE_GENOME ( SAMTOOLS_FAIDX.out.fai )

    emit:
    dot_genome      = GENERATE_GENOME.out.dotgenome

    versions        = ch_versions.ifEmpty(null)

}