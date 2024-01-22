process CHUNKFASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::pyfasta=0.5.2-1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/pyfasta:0.5.2--py_1' :
    'biocontainers/pyfasta:0.5.2--py_1' }"

    input:
    tuple val(meta), path('input.fasta')
    val(number_of_chunks)

    output:
    tuple val(meta), path('*.fasta'), emit: fasta
    path "versions.yml"             , emit: versions

    script:
    def VERSION = '0.5.2' // Tool does not report version
    // This should be abstracted outside of the container to
    // stop it spinning up in the first place,
    // however dsl2 can't do comparisons with channels which makes it harder
    """
    if [ $number_of_chunks -le 1 ]; then
        mv input.fasta ${meta.id}_whole.fasta
    else
        pyfasta split -n $number_of_chunks input.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyfasta: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '0.5.2' // Tool does not report version
    """
    touch ${meta.id}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pyfasta: $VERSION
    END_VERSIONS
    """
}
