process GRAPHOVERALLCOVERAGE {
    tag "$meta.id"
    label "process_single"

    container 'quay.io/sanger-tol/cramfilter_bwamem2_minimap2_samtools_perl:0.001-c1'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "GRAPHOVERALLCOVERAGE module does not support Conda. Please use Docker / Singularity instead."
    }

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.part") , emit: part
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    graph_overall_coverage.pl $bed > ${prefix}.part

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        graph_overall_coverage.pl: \$(graph_overall_coverage.pl --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.part

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        graph_overall_coverage.pl: \$(graph_overall_coverage.pl --version)
    END_VERSIONS
    """
}
