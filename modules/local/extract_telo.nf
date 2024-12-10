process EXTRACT_TELO {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*bed" )   , emit: bed
    tuple val( meta ), file("*bedgraph"), emit: bedgraph
    path "versions.yml"                 , emit: versions

    shell:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    $/
    cat "${file}" |awk '{print $2"\t"$4"\t"$5}'|sed 's/>//g' > ${prefix}_telomere.bed
    cat "${file}" |awk '{print $2"\t"$4"\t"$5"\t"((($5-$4)<0)?-($5-$4):($5-$4))}' | sed 's/>//g' > ${prefix}_telomere.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    /$

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_telomere.bed
    touch ${prefix}_telomere.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
