process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/sanger-tol/pretext:0.0.8-yy5-c1"

    input:
    tuple val(meta),    path(pretext_file)
    tuple val(gap),     path(gap_file)
    tuple val(cov),     path(coverage)
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
    if [[ -s "${coverage}" ]];
    then
        echo "PROCESSING COVERAGE..."
        bigWigToBedGraph ${coverage} /dev/stdout | PretextGraph ${args} -i ${pretext_file} -n "coverage" -o coverage.pretext.part
    else
        echo "SKIPPING COVERAGE"
        mv ${pretext_file} coverage.pretext.part
    fi

    if [[ -s "${repeat_density}" ]];
    then
        echo "PROCESSING REPEAT_DENSITY..."
        bigWigToBedGraph  ${repeat_density} /dev/stdout | PretextGraph ${args} -i coverage.pretext.part -n "repeat_density" -o repeat.pretext.part
    else
        echo "SKIPPING REPEAT_DENSITY"
        mv coverage.pretext.part repeat.pretext.part
    fi

    if [[ -s "${gap_file}" ]]; then
        echo "Processing GAP file..."
        cat "${gap_file}" | PretextGraph ${args} -i repeat.pretext.part -n "${gap.ft}" -o gap.pretext.part
        input_file="gap.pretext.part"
    else
        input_file="repeat.pretext.part"
    fi

    if [[ -s "${telomere_file}" ]]; then
        echo "Processing TELO file..."
        cat "${telomere_file}" | PretextGraph ${args} -i "${input_file}" -n "${telo.ft}" -o "${prefix}.pretext"
    else
        mv "${input_file}" "${prefix}.pretext"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext.* Version //;')
        PretextMap: \$(PretextMap | grep "Version" | sed 's/Pretext.* Version//;')
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
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext* Version //;')
        PretextMap: \$(PretextMap | grep "Version" | sed 's/PretextMap Version//;')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """
}
