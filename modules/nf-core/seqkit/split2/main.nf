process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("**/*.{fa,fasta,gz}"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    if ( !meta.single_end && meta.file_type == "fastq" ) {

        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            --read1 ${reads[0]} \\
            --read2 ${reads[1]} \\
            --out-dir ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            $reads \\
            --out-dir ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( !meta.single_end && meta.file_type == "fastq" ) {
        """
        mkdir -p ${prefix}
        echo "" | gzip > ${prefix}/${reads[0]}
        echo "" | gzip > ${prefix}/${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        mkdir -p ${prefix}
        echo "" | gzip > ${prefix}/${reads[0]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
