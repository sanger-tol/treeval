process CONCATMUMMER {
<<<<<<< HEAD
    tag "${meta.id}.mummer"
    label "process_medium"

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the mummer2bed process. Please use docker or singularity containers."
    }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
            'ubuntu:20.04' }"
=======
    tag "${meta} -> ${meta.id}.mummer"
    label "process_small"
>>>>>>> 24a3749 (All modules added)

    input:
    tuple val(meta), path(coords)

    output:
    tuple val(meta), path("*.mummer"), emit: mummer
<<<<<<< HEAD
    path "versions.yml", emit: versions
=======
>>>>>>> 24a3749 (All modules added)

    script:
    """
    cat $coords > ${meta.id}.mummer
<<<<<<< HEAD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ubuntu: \$(ubuntu --version | sed 's/Ubuntu //g')
=======
>>>>>>> 24a3749 (All modules added)
    """
}
