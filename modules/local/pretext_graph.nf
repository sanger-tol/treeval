process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextgraph=0.0.6 bioconda::ucsc-bigwigtobedgraph=448"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-077b852b8b5440d395ad23f9f24f50c943390a84:da499c75fec554e81f4847c4fa8b6b167afbe3bf-0':
        'biocontainers/mulled-v2-077b852b8b5440d395ad23f9f24f50c943390a84:da499c75fec554e81f4847c4fa8b6b167afbe3bf-0' }"

    input:
    tuple val(meta) , path(pretext_file)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(repeat_density)
    tuple val(meta4), path(log_coverage)
    tuple val(meta5), path(gap_file)
    tuple val(meta6), path(telomere_file)

    output:
    tuple val(meta), path("*.pretext")  , emit: pretext
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    bigWigToBedGraph ${coverage} /dev/stdout | PretextGraph -i ${pretext_file} -n "coverage" -o coverage.pretext.part

    bigWigToBedGraph  ${repeat_density} /dev/stdout | PretextGraph -i coverage.pretext.part -n "repeat_density" -o repeat.pretext.part

    bigWigToBedGraph  ${log_coverage} /dev/stdout | PretextGraph -i repeat.pretext.part -n "log2_coverage" -o log.pretext.part

    if [[ ${meta5.sz} -ge 1 && ${meta6.sz} -ge 1 ]]
    then
        echo "GAP AND TELO have contents!"
        cat ${gap_file} | PretextGraph -i log.pretext.part -n "${meta5.ft}" -o gap.pretext.part
        cat ${telomere_file} | awk -v OFS='\t' {\$4 *= 1000; print} | PretextGraph -i gap.pretext.part -n "${meta6.ft}" -o ${prefix}.pretext

    elif [[ ${meta5.sz} -ge 1 && ${meta6.sz} < 0 ]]
    then
        echo "GAP file has contents!"
        cat ${gap_file} | PretextGraph -i ${prefix}.pretext.part -n "${meta5.ft}" -o ${prefix}.pretext

    elif [[ ${meta5.sz} < 0 && ${meta6.sz} -ge 1 ]]
    then
        echo "TELO file has contents!"
        cat ${telomere_file} | awk -v OFS='\t' {\$4 *= 1000; print} | PretextGraph -i ${pretext_file} -n "${meta6.ft}" -o ${prefix}.pretext

    else
        echo "NO GAP OR TELO FILE WITH CONTENTS - renaming part file"
        mv log.pretext.part ${prefix}.pretext
    fi

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
    touch ${prefix}.pretext

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextmap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version //g')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """
}
