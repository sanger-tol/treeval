include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { GENERATE_GENOME_FILE  } from '../../modules/local/generate_genome_file'

workflow GENERATE_GENOME {
    take:
    assembly_id
    reference_file

    main:
    ch_versions     = Channel.empty()

<<<<<<< HEAD
    reference_file
        .combine( assembly_id )
        .map { it ->
            tuple ([id: it[1]],
                    it[0])
        }
=======
    TO_FILE ( assembly_id, reference_file)

    TO_FILE.out.file_path
        .map( it -> 
                tuple(
                    [id: it[0]],
                    it[1]
                )
        )
>>>>>>> 387355f (Adding GENERATE_GENOME subworkflow to main)
        .set { to_samtools }

    SAMTOOLS_FAIDX ( to_samtools )
    ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    GENERATE_GENOME_FILE ( SAMTOOLS_FAIDX.out.fai )
 
    emit:
    dot_genome      = GENERATE_GENOME_FILE.out.dotgenome
    reference_tuple = to_samtools

    versions        = ch_versions.ifEmpty(null)
}
