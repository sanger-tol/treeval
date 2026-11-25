process PRETEXT_GRAPH {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::pretextgraph=0.1.0"
    container "quay.io/biocontainers/pretextgraph:0.1.0--hc52dbad_0"

    input:
    tuple val(meta), path(pretext_file), path(tracks, arity: '0..*'), val(track_names)

    output:
    tuple val(meta), path("*.pretext"), emit: pretext
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def track_list = tracks instanceof List ? tracks : (tracks ? [tracks] : [])
    def name_list  = track_names instanceof List ? track_names : (track_names != null ? [track_names] : [])

    def track_checks = track_list.withIndex().collect { track_file, i ->
        def track_name = (i < name_list.size() && name_list[i]) ? name_list[i] : track_file.baseName
        """
        if [ -s "${track_file}" ]; then
            echo "Adding track ${track_name} from ${track_file}"
            track_args="\${track_args} -t\${idx} ${track_file} -n\${idx} ${track_name}"
            idx=\$((idx + 1))
        else
            echo "Skipping track ${track_name} (empty or missing)"
        fi
        """
    }.join('\n')

    """
    track_args=""
    idx=1
    ${track_checks}

    if [ -n "\$track_args" ]; then
        PretextGraph \\
            ${args} \\
            -i ${pretext_file} \\
            \$track_args \\
            -o ${prefix}.pretext
    else
        echo "No accessory tracks to ingest; copying input pretext"
        cp ${pretext_file} ${prefix}.pretext
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext.* Version //;')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pretext
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextGraph: \$(PretextGraph | grep "Version" | sed 's/Pretext.* Version //;')
    END_VERSIONS
    """
}
