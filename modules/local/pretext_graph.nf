process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::pretextgraph=0.0.6 bioconda::ucsc-bigwigtobedgraph=448"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-077b852b8b5440d395ad23f9f24f50c943390a84:da499c75fec554e81f4847c4fa8b6b167afbe3bf-0':
        'biocontainers/mulled-v2-077b852b8b5440d395ad23f9f24f50c943390a84:da499c75fec554e81f4847c4fa8b6b167afbe3bf-0' }"

    input:
    tuple val(meta),    path(pretext_file)
    tuple val(gap),     path(gap_file)
    tuple val(cov),     path(coverage)
    tuple val(log),     path(log_coverage)
    tuple val(telo),    path(telomere_file)
    tuple val(rep),     path(repeat_density)

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

    bigWigToBedGraph  ${log_coverage} /dev/stdout | awk -v OFS="\t" '{ if (\$4 < 4) {\$4 *= 1000} else {\$4 *= 100} ; print}' | PretextGraph -i repeat.pretext.part -n "log_coverage" -o log.pretext.part

    if [[ ${gap.sz} -ge 1 && ${telo.sz} -ge 1 ]]
    then
        echo "GAP AND TELO have contents!"
        cat ${gap_file} | PretextGraph -i log.pretext.part -n "${gap.ft}" -o gap.pretext.part
        cat ${telomere_file} | awk -v OFS='\t' '{\$4 *= 1000; print}' | PretextGraph -i gap.pretext.part -n "${telo.ft}" -o ${prefix}.pretext

    elif [[ ${gap.sz} -ge 1 && ${telo.sz} -eq 0 ]]
    then
        echo "GAP file has contents!"
        cat ${gap_file} | PretextGraph -i log.pretext.part -n "${gap.ft}" -o ${prefix}.pretext

    elif [[ ${gap.sz} -eq 0 && ${telo.sz} -ge 1 ]]
    then
        echo "TELO file has contents!"
        cat ${telomere_file} | awk -v OFS='\t' '{\$4 *= 1000; print}' | PretextGraph -i log.pretext.part -n "${telo.ft}" -o ${prefix}.pretext

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
