//
// Check for synteny by aligning to fasta to reference genomes.
//
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'

workflow SYNTENY {
    take:
        reference_tuple
        synteny_path
        assembly_classT

    main:
    ch_versions = Channel.empty()

    // If no reference genomes returned for class then no .PAF outputted.
    //ch_genomes = Channel.fromPath("${synteny_path}/${assembly_classT}/*.fasta")
    
    process GET_SYNTENY_GENOMES {
        tag "${sample_id}"
        label "process_small"

        input:
        path ( synteny_path )
        val ( assembly_classT )
        
        output:
        path ( '*fasta' ), emit: genome_path
        
        shell:
        def synteny_filepath = synteny_path.toString()
        """
        if [ -f "${synteny_path}/${assembly_classT}/*.fasta" ] ; then
            cp "$synteny_filepath/${assembly_classT}/*.fasta" "./*.fasta"
        else
            cat > empty.fasta
        fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GET_SYNTENY_GENOMES: BASH
        END_VERSIONS
        """
    }

    GET_SYNTENY_GENOMES(synteny_path, assembly_classT)

    ch_genome_paths = GET_SYNTENY_GENOMES.out.genome_path
    .branch {
        genome_paths -> 
        empty : genome_paths == "./empty.fasta"
            return genome_paths
        refs : genome_paths != "./empty.fasta"
            return genome_paths
    }
    ch_genome_paths.empty.view()
    ch_genome_paths.refs.view()

    MINIMAP2_ALIGN(reference_tuple, ch_genome_paths.refs, false, true, false)

    ch_paf = MINIMAP2_ALIGN.out.paf
    ch_versions = MINIMAP2_ALIGN.out.versions
    
    emit:
    ch_paf
    versions = ch_versions
}