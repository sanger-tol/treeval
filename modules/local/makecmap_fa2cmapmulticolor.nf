process MAKECMAP_FA2CMAPMULTICOLOR {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::perl=5.26.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl:5.26.2' :
        'quay.io/biocontainers/perl:5.26.2' }"

    input:
    tuple val(meta), path(fasta)
    val enzyme

    output:
    tuple val(meta), path("*.cmap"), emit: cmap
    path("*key.txt")               , emit: cmapkey
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fa2cmap_multi_color.pl -i $fasta -e $enzyme 1 $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        fa2cmap_multi_color.pl: \$(fa2cmap_multi_color.pl -v)
    END_VERSIONS
    """
}
