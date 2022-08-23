process CHUNKFASTA {
    tag "$chunkfasta"

    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
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
