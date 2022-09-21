process CHUNKFASTA {
    tag "${meta.id}"

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the chunkfasta process. Please use docker or singularity containers."
    }

    tag "$chunkfasta"

    container "quay.io/biocontainers/pyfasta:0.5.2--py_1"

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
