include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { GENERATE_GENOME_FILE  } from '../../modules/local/genome_file_generator'

workflow GENERATE_GENOME {
    main:
    ch_versions     = Channel.empty()

    SAMTOOLS_FAIDX ( [[params.assembly.sample], params.reference] )
    ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GENERATE_GENOME_FILE ( SAMTOOLS_FAIDX.out.fai )

    emit:
    dot_genome      = GENERATE_GENOME_FILE.out.dotgenome

    versions        = ch_versions.ifEmpty(null)

}