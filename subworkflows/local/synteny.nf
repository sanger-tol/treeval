//
// Check for synteny by aligning to fasta to reference genomes.
//
include { MINIMAP2_ALIGN        } from '../../modules/nf-core/modules/minimap2/align/main'
include { GET_SYNTENY_GENOMES   } from '../../modules/local/get_synteny_genomes'

workflow SYNTENY {
    take:
        reference_tuple
        synteny_path
        assembly_classT

    main:
    ch_versions = Channel.empty()

    GET_SYNTENY_GENOMES(synteny_path, assembly_classT)

    GET_SYNTENY_GENOMES.out.genome_path
        .flatten()
        .branch { data -> 
            run: !data.toString().contains("empty")
            skip: data.toString().contains("empty")
        }
        .set { mm_intermediary }

    reference_tuple
        .combine(mm_intermediary.run)
        .map { meta, fa, ref ->
            tuple(meta, fa, ref, false, true, false)
            }
        .set { mm_input }

    MINIMAP2_ALIGN( mm_input.map { [it[0], it[1]] },
                    mm_input.map { it[2] },
                    mm_input.map { it[3] },
                    mm_input.map { it[4] },
                    mm_input.map { it[5] } )

    ch_paf = MINIMAP2_ALIGN.out.paf
    ch_versions = MINIMAP2_ALIGN.out.versions
    
    emit:
    ch_paf
    versions = ch_versions
}