process SEQKIT_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path('input.fasta')
    val(number_of_chunks)

    output:
    tuple val(meta), path('*.{fa,fasta}'), emit: fasta
    path "versions.yml",            emit: versions

    script:
    // This should be abstracted outside of the container to
    // stop it spinning up in the first place,
    // however dsl2 can't do comparisons with channels which makes it harder
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ $number_of_chunks -le 1 ]; then
        mv input.fasta ${prefix}_whole.fa
    else
        seqkit split input.fasta -p $number_of_chunks --by-part-prefix ${prefix}_ -O ./
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}
