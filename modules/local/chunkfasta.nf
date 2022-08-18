process CHUNKFASTA {
    tag "$chunkfasta"

    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "quay.io/biocontainers/pyfasta:0.5.2--py_1"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.fa'), emit: fas
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/treeval/bin/
    //def number_of_chunks = fasta.size()
    """
    file_size = os.path.getsize(${fasta})
    file_size_in_Gb = filesize/1073741824
    number_of_chunks = ceil(file_size_in_Gb)
    
    pyfasta split -n 6 $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
