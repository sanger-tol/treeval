process MERGE_BAM {
    tag " ${meta.id} "
    label "process_medium"

    container ""

    script:
    """
    samtools merge merged.bam *.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}": 
        MERGE_BAM: $version
    END_VERSIONS
    """

}