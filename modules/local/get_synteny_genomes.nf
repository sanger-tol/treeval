    process GET_SYNTENY_GENOMES {
        tag "${assembly_classT}"
        label "process_low"

        conda (params.enable_conda ? "conda-forge::coreutils=9.1" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

        input:
        val ( synteny_path )
        val ( assembly_classT )
        
        output:
        path ( '*fasta' ), emit: genome_path
        
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
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            GET_SYNTENY_GENOMES: BASH
        END_VERSIONS
        """
    }