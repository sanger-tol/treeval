
process PRETEXTMAP {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextmap=0.1.9 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1' :
        'biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(input)
    path fasta

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION         = "0.1.9"
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def reference       = fasta ? "--reference ${fasta}" : ""
    def pretext_path    = "${projectDir}/bin/PretextMap/bin/PretextMap"

    """
    if [[ $input == *.pairs.gz ]]; then
        zcat $input | ${pretext_path} \\
            $args \\
            -o ${prefix}.pretext
    else
        samtools \\
            view \\
            $reference \\
            -h \\
            $input | ${pretext_path} \\
            $args \\
            -o ${prefix}.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    def VERSION         = "0.1.9"
    def prefix          = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: $VERSION
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
