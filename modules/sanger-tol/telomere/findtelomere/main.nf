process TELOMERE_FINDTELOMERE {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/sanger-tol/telomere:0.0.4-c2'

    input:
    tuple val(meta), path(reference), val(telomereseq)

    output:
    tuple val(meta), path("${prefix}.telomere"), emit: telomere
    // FindTelomereWindows always writes combined motif BED; strand BEDs only with --split.
    tuple val(meta), path("${prefix}.telomere.bed"), emit: telomere_bed
    tuple val(meta), path("${prefix}.fwd.telomere.bed"), emit: telomere_bed_fwd, optional: true
    tuple val(meta), path("${prefix}.rev.telomere.bed"), emit: telomere_bed_rev, optional: true
    // Combined density windows always written as `prefix.windows`; strand windows only with --split.
    tuple val(meta), path("${prefix}.windows"), emit: windows_all
    tuple val(meta), path("${prefix}.fwd.windows"), emit: windows_fwd, optional: true
    tuple val(meta), path("${prefix}.rev.windows"), emit: windows_rev, optional: true
    tuple val("${task.process}"), val('java'), eval("java -version 2>&1 | head -n 1 | cut -d '\"' -f2"), topic: versions, emit: versions_java
    // find_telomere has no --version; pin to container tag (bump when `container` changes)
    tuple val("${task.process}"), val('find_telomere'), val('0.0.4'), topic: versions, emit: versions_find_telomere

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def max_heap_size_mega = (task.memory.toMega() * 0.9).intValue()
    def max_stack_size_mega = 999 //most java jdks will not allow Xss > 1GB, so fixing this to the allowed max

    """
    find_telomere ${args} ${reference} ${telomereseq} | awk '{print \$1"\\t"\$(NF-4)"\\t"\$(NF-3)"\\t"\$(NF-2)"\\t"\$(NF-1)"\\t"\$NF}' - > ${prefix}.telomere

    java \\
        -Xmx${max_heap_size_mega}M \\
        -Xss${max_stack_size_mega}M  \\
        -cp /opt/telomere/telomere.jar \\
        FindTelomereWindows \\
        ${prefix}.telomere \\
        ${args2}
    """

    stub:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FINDTELOMERE module does not support Conda. Please use Docker / Singularity instead."
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def args2 = task.ext.args2 ?: ''

    """
    printf "stub\\n" > ${prefix}.telomere
    printf "stub\\n" > ${prefix}.telomere.bed
    printf "stub\\n" > ${prefix}.windows
    if [[ "${args2}" == *--split* ]]; then
        printf "stub\\n" > ${prefix}.fwd.telomere.bed
        printf "stub\\n" > ${prefix}.rev.telomere.bed
        printf "stub\\n" > ${prefix}.fwd.windows
        printf "stub\\n" > ${prefix}.rev.windows
    fi
    """

}
