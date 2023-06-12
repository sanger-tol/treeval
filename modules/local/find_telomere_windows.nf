process FIND_TELOMERE_WINDOWS {
    tag "${meta.id}"
    label "process_low"

    conda "bioconda::java-jdk=8.0.112"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val( meta ), path( file )

    output:
    tuple val( meta ), file( "*.windows" ) , emit: windows
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def telomere_jar = task.ext.telomere_jar ?: ''
    def telomere_jvm_params = task.ext.telomere_jvm_params ?: ''
    """
    java ${telomere_jvm_params} -cp ${projectDir}/bin/${telomere_jar} FindTelomereWindows $file 99.9 > ${prefix}.windows

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomere: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def telomere = task.ext.telomere ?: ''
    """
    touch ${prefix}.windows

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomere: $VERSION
    END_VERSIONS
    """

}
