process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextgraph=0.6.3 bioconda::ucsc-bedtobigbed=377"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextgraph%3A0.0.6--h9f5acd7_2':
        'biocontainers/pretextgraph:0.0.6' }"

    input:
    tuple val(meta) , path(pretext_file)
    tuple val(meta2), path(gap)
    tuple val(meta3), path(coverage)
    tuple val(meta4), path(telomere)
    tuple val(meta5), path(repeat_density)

    output:
    tuple val(meta), path("*_updated.pretext")  , emit: pretext
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bigWigToBedGraph ${coverage} /dev/stdout | PretextGraph -i ${pretext_file} -n "coverage"

    bigWigToBedGraph ${repeat_density} /dev/stdout | PretextGraph -i ${pretext_file} -n "repeat_density"

    cat ${telomere} | awk -v OFS="\t" '{\$4 *= 1000; print}' | PretextGraph -i ${pretext_file} -n "telomere"

    cat ${gap} | PretextGraph -i ${pretext_file} -n "gap"

    mv ${pretext_file} ${prefix}_updated.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextgraph: \$(PretextGraph | grep "Version" | sed 's/PretextGraph Version //g')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '377' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_updated.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version //g')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """
}
