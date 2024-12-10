process GENERATE_CRAM_CSV {
    tag "${meta.id}"
    label 'process_tiny'

    container 'quay.io/sanger-tol/cramfilter_bwamem2_minimap2_samtools_perl:0.001-c1'

    input:
    tuple val(meta), path(crampath)


    output:
    tuple val(meta), path('*.csv'), emit: csv
    path "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_cram_csv.sh $crampath ${prefix}_cram.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_cram_csv: \$(generate_cram_csv.sh -v)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generate_cram_csv: \$(generate_cram_csv.sh -v)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
