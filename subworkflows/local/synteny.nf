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
        .set { mm_input }

    mm_input.run.view()

    MINIMAP2_ALIGN(reference_tuple, mm_input.run, false, true, false)

    ch_paf = MINIMAP2_ALIGN.out.paf
    ch_versions = MINIMAP2_ALIGN.out.versions
    
    emit:
    ch_paf
    versions = ch_versions
}