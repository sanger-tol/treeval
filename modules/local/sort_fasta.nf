process SORT_FASTA {
    tag "$meta.id"
    label "process_single"

    conda "conda-forge::ruby=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ruby:2.2.3--1' :
        'biocontainers/ruby:2.2.3--1' }"

    input:
    tuple val(meta), path(input_assembly)
    val (output_prefix)
    val (scaffold_length_cutoff)

    output:
    tuple val(meta), path("*.fa"), emit: fa
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "2.2.3"
    """
    sort_fasta.rb -f ${input_assembly} -o ${output_prefix}.MicroFinder.order.tsv -l ${scaffold_length_cutoff} > ${output_prefix}.MicroFinder.ordered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ruby: $VERSION
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "2.2.3"
    """
    touch ${output_prefix}.MicroFinder.ordered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ruby: $VERSION
    END_VERSIONS
    """
}
