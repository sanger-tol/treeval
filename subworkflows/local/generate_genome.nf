include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { GENERATE_GENOME_FILE  } from '../../modules/local/generate_genome_file'

workflow GENERATE_GENOME {
    take:
    assembly_id
    reference_file

    main:
    ch_versions     = Channel.empty()

    reference_file
        .combine( assembly_id )
        .map { it ->
            tuple ([id: it[1]],
                    it[0])
        }
        .set { to_samtools }

    SAMTOOLS_FAIDX ( to_samtools )
    ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GENERATE_GENOME_FILE ( SAMTOOLS_FAIDX.out.fai )
 
    emit:
    dot_genome      = GENERATE_GENOME_FILE.out.dotgenome
    reference_tuple = to_samtools

    versions        = ch_versions.ifEmpty(null)
}

