process CHECK_RG {
    tag "${meta.id}"
    label 'process_tiny'

    container 'quay.io/sanger-tol/cramfilter_bwamem2_minimap2_samtools_perl:0.001-c1'

    input:
    tuple val(meta), path(crampath)
    path outdir
 
    output:
    tuple val(meta), path("*.cram")  , emit: newcram
    tuple val(meta), path("*.crai")  , emit: newcrai
    path "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    split_rg_cram.sh $crampath $outdir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cram
    touch ${meta.id}.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
