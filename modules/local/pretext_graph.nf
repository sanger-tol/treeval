process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextgraph=0.0.6 bioconda::ucsc-bigwigtobedgraph=448"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4a53f3db479fca2f9afe0c7dc59a6fe8a053b096:298ec4f1a4dd566f5ec43ff2fd96b8502f9613fd-0':
        'biocontainers/mulled-v2-4a53f3db479fca2f9afe0c7dc59a6fe8a053b096:298ec4f1a4dd566f5ec43ff2fd96b8502f9613fd-0' }"

    input:
    tuple val(meta) , path(pretext_file)
    tuple val(meta2), path(gap)
    tuple val(meta3), path(telomere)
    tuple val(meta4), path(coverage)
    tuple val(meta5), path(repeat_density)

    output:
    tuple val(meta), path("*_updated.pretext")  , emit: pretext
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bigWigToBedGraph ${coverage} /dev/stdout | PretextGraph -i ${pretext_file} -n "coverage"

    bigWigToBedGraph ${repeat_density} /dev/stdout | PretextGraph -i *.pretext -n "repeat_density"

    cat ${telomere} | awk -v OFS="\t" '{\$4 *= 1000; print}' | PretextGraph -i *.pretext -n "telomere"

    cat ${gap} | PretextGraph -i *.pretext -n "gap"

    mv *.pretext ${prefix}_updated.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextgraph: \$(PretextGraph | grep "Version" | sed 's/PretextGraph Version //g')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '448' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_updated.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version //g')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """
}
