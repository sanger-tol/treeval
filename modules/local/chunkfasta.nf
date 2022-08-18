process CHUNKFASTA {
<<<<<<< HEAD
    tag "${meta.id}"

    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used when using the chunkfasta process. Please use docker or singularity containers."
    }
=======
    tag "$chunkfasta"

    //conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
>>>>>>> 46fc808 (Start selfcomp workflow)
    container "quay.io/biocontainers/pyfasta:0.5.2--py_1"

    input:
    tuple val(meta), path(fasta)
<<<<<<< HEAD
    val(number_of_chunks)
=======
>>>>>>> 46fc808 (Start selfcomp workflow)

    output:
    tuple val(meta), path('*.fa'), emit: fas
    path "versions.yml", emit: versions

<<<<<<< HEAD
    script:
    """
    pyfasta split -n $number_of_chunks $fasta
=======
    script: // This script is bundled with the pipeline, in nf-core/treeval/bin/
    //def number_of_chunks = fasta.size()
    """
    file_size = os.path.getsize(${fasta})
    file_size_in_Gb = filesize/1073741824
    number_of_chunks = ceil(file_size_in_Gb)
    
    pyfasta split -n 6 $fasta
>>>>>>> 46fc808 (Start selfcomp workflow)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
