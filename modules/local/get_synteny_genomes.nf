process GET_SYNTENY_GENOMES {
    tag "${defined_class}"
    label "process_single"

    conda "conda-forge::coreutils=9.1"

    input:
    val ( synteny_path )
    val ( defined_class )

    output:
    path ( '*fasta' )   , emit: genome_path
    path "versions.yml" , emit: versions

    shell:
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    '''
    VERSION="9.1"

    if [[ -d "!{synteny_path}!{defined_class}/" ]]
    then
        if [[ "$(ls -A !{synteny_path}!{defined_class}/)" ]]
        then
            echo "Found Some Sytenic Genomes"
            cp !{synteny_path}!{defined_class}/*.fasta .
        else
            echo "Directory exists but is empty"
            touch empty.fasta
        fi
    else
        echo "Directory doesn't exist"
        touch empty.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
        bash: $(echo $(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    '''

    stub:
    def VERSION     = "9.1" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch empty.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version | grep -Eo 'version [[:alnum:].]+' | sed 's/version //'))
        coreutils: $VERSION
    END_VERSIONS
    """
}
