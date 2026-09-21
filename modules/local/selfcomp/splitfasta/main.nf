process SELFCOMP_SPLITFASTA {
    tag "${meta.id}"
    label "process_single"

    conda "conda-forge::perl-bioperl=1.7.8-1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.8--hdfd78af_1' :
        'biocontainers/perl-bioperl:1.7.8--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fa
    path("*.agp")                , emit: agp
    tuple val("${task.process}"), val('split_genomes_for_ensembl.pl'), eval("split_genomes_for_ensembl.pl --version"), topic: versions, emit: versions_split_genomes_for_ensembl
    tuple val("${task.process}"), val('perl'), eval("perl --version | sed -n 's/.*(v\\([0-9.]\\+\\)).*/\\1/p'"), topic: versions, emit: versions_perl
    tuple val("${task.process}"), val('perl-bioperl'), val("1.7.8-1"), topic: versions, emit: versions_perlbioperl

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "1.7.8-1"
    """
    split_genomes_for_ensembl.pl $fasta ${prefix}_windowed.fa ${prefix}_split.agp
    """

    stub:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def VERSION     = "1.7.8-1"
    """
    touch ${prefix}_split.agp
    touch ${prefix}_split.fa
    """
}
