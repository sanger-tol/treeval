process FILTER_BLAST {
    tag "${meta.id} - ${meta.type}"
    label "process_medium"

    def version = 'v1.0.0'

    /*     conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/python:3.8.3' :
            'quay.io/biocontainers/python:3.8.3' }" */

    input:
    tuple val( meta ), file( concat_blast_out )

    output:
    tuple val( meta ), file( "${meta.id}-${meta.type}-*.tsv")   , emit: final_tsv
    path "versions.yml"                                         , emit: versions

    script:
    def filt_percent = task.ext.args ?: 90.00
    // HOW DO I INSTALL PANDAS INTO THIS? MAKE OWN CONTAINER?
    """
    /software/grit/conda/envs/Damon_project/bin/python3 $projectDir/bin/filter_blast.py $meta.id $meta.type $concat_blast_out $filt_percent

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blast: $version
    END_VERSIONS
    """
}
