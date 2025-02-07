process MICROFINDER_FILTER {
    tag "$meta.id"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'docker.io/ubuntu:20.04' }"

    input:
    tuple val(meta), path(microfinder_gff)
    tuple val(meta), path(input_assembly)
    val (output_prefix)
    val (scaffold_length_cutoff)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    awk '${scaffold_length_cutoff}== "mRNA"' ${microfinder_gff} | grep -w "Rank=1" | cut -f1,9 | tr ";" "\t" | cut -f1,4,5,6,7 | sed 's/Identity=//g' | sed 's/Positive=//g' | awk '${output_prefix} >= 0.7' | cut -f1 | sort | uniq -c | sort -k1,1nr | awk '{print ${output_prefix} "\t" ${input_assembly}}' > ${output_prefix}.MicroFinder.order.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${output_prefix}.MicroFinder.order.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: $VERSION
    END_VERSIONS
    """
}
