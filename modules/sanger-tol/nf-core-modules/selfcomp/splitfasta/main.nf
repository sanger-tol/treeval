process SELFCOMP_SPLITFASTA {
    tag "$meta.id"
    label 'process_medium'

    def version = '0.001-c3'

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the splitfasta process. Please use docker or singularity containers."
    }
    container "quay.io/sanger-tol/splitfasta:${version}"

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
        selfcomp_splitfasta: ${version}
    END_VERSIONS
    """
}

