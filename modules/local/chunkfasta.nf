process CHUNKFASTA {
    tag "${meta.id}"
    label "process_medium"

    conda "conda-forge::pyfasta=0.5.2-1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pyfasta:0.5.2--py_1' :
    'quay.io/biocontainers/pyfasta:0.5.2--py_1' }"

    input:
    tuple val(meta), path(fasta)
    val(number_of_chunks)

    output:
    tuple val(meta), path('*.fa'), emit: fas
    path "versions.yml", emit: versions

    script:
    """
    pyfasta split -n $number_of_chunks $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
