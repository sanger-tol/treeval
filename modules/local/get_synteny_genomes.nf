process GET_SYNTENY_GENOMES {
    tag "${assembly_classT}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
    'ubuntu:20.04' }"

    input:
    val ( synteny_path )
    val ( assembly_classT )

    output:
    path ( '*fasta' )   , emit: genome_path
    path "versions.yml" , emit: versions

    script:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ ! -d ${synteny_path}${assembly_classT}/ ] || [ -z "\$(ls -A ${synteny_path}${assembly_classT}/)" ]
    then
        echo "Directory is empty or doesn't exist"
        touch empty.fasta
    else
        cp ${synteny_path}${assembly_classT}/*.fasta .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch empty.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    """
}
