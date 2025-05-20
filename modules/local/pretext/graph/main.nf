process PRETEXT_GRAPH {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/sanger-tol/pretext:0.0.9-yy5-c2"

    input:
    tuple val(meta),        path(pretext_file)
    path(gap_file,          stageAs: 'gap_file.bed')
    path(coverage,          stageAs: 'coverage.bw')
    path(telomere_file,     stageAs: 'telomere.bed')
    path(repeat_density,    stageAs: 'repeat_density.bw')

    output:
    tuple val(meta), path("*.pretext")  , emit: pretext
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT GRAPH module does not _currently_ support Conda. Please use Docker / Singularity instead."
    }

    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def UCSC_VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // Using single [ ] as nextflow will use sh where possible not bash
    """

    echo "PROCESSING ESSENTIAL FILES"

    if [ -s "${coverage}" ]; then
        echo "PROCESSING COVERAGE..."
        bigWigToBedGraph ${coverage} /dev/stdout | PretextGraph ${args} -i ${pretext_file} -n "coverage" -o coverage.pretext.part
    else
        echo "SKIPPING COVERAGE"
        mv ${pretext_file} coverage.pretext.part
    fi

    if [ -s "${repeat_density}" ]; then
        echo "PROCESSING REPEAT_DENSITY..."
        bigWigToBedGraph  ${repeat_density} /dev/stdout | PretextGraph ${args} -i coverage.pretext.part -n "repeat_density" -o repeat.pretext.part
    else
        echo "SKIPPING REPEAT_DENSITY"
        mv coverage.pretext.part repeat.pretext.part
    fi

    echo "NOW PROCESSING NON-ESSENTIAL files"

    input_file="repeat.pretext.part"

    if [ -s "${gap_file}" ]; then
        echo "Processing GAP file..."
        cat "${gap_file}" | PretextGraph ${args} -i repeat.pretext.part -n "gap" -o gap.pretext.part
        input_file="gap.pretext.part"
    fi

    if [ -s "${telomere_file}" ]; then
        echo "Processing TELO file..."
        cat "${telomere_file}" | PretextGraph ${args} -i "\$input_file" -n "telomere" -o "${prefix}.pretext"
    else
        mv "\$input_file" "${prefix}.pretext"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext.* Version //;')
        PretextMap: \$(PretextMap | grep "Version" | sed 's/Pretext.* Version//;')
        bigWigToBedGraph: ${UCSC_VERSION}
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PRETEXT GRAPH module does not _currently_ support Conda. Please use Docker / Singularity instead."
    }

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
