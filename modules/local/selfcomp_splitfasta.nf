process SELFCOMP_SPLITFASTA {
    tag "$meta.id"
    label "process_medium"

    conda "conda-forge::perl-bioperl=1.7.8-1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.8--hdfd78af_1' :
        'quay.io/biocontainers/perl-bioperl:1.7.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fa
    path("*.agp")                , emit: agp
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_genomes_for_ensembl.pl $fasta ${prefix}_split.fa ${prefix}_split.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        perl-bioperl: 1.7.8-1
        split_genomes_for_ensembl.pl: \$(split_genomes_for_ensembl.pl --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_split.agp
    touch ${prefix}_split.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/^.*perl //; s/Using.*\$//')
        perl-bioperl: 1.7.8-1
        split_genomes_for_ensembl.pl: \$(split_genomes_for_ensembl.pl --version)
    END_VERSIONS
    """
}

