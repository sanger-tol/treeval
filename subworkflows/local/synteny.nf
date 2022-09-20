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
        val ( synteny_path )
        val ( assembly_classT )
        
        output:
        path "*.fasta", emit: file_path
        
        script:
        """
        if [ ! -d ${synteny_path}/${assembly_classT}/ ] || [ -z "\$(ls -A ${synteny_path}/${assembly_classT}/)" ]
        then
            echo "Directory is empty or doesn't exist"
            touch empty.fasta
        else
            for i in "${synteny_path}/${assembly_classT}/*.fasta"
                do
                    cp \$i .
                done
        fi
        """
    }

    GET_SYNTENY_GENOMES(synteny_path, assembly_classT)

    GET_SYNTENY_GENOMES.out.file_path.flatten()
        .branch { data -> 
            run: !data.contains("empty")
            skip: data.contains("empty")
        }
        .set { ch_genome_paths }

    MINIMAP2_ALIGN(reference_tuple, ch_genome_paths.run, false, true, false)

    ch_paf = MINIMAP2_ALIGN.out.paf
    ch_versions = MINIMAP2_ALIGN.out.versions
    
    emit:
    ch_paf
    versions = ch_versions
}