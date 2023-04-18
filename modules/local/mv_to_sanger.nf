process MV_TO_SANGER {
    input:
    val( ref_name )                         // iyTipFemo1
    val( dbVersion )                        // 1
    tuple val( meta ), path( punchlist )    // [[ meta.id ], punchlist.bed ]

    output:
    tuple   val( meta ), file( punchlist )

    script:
    """
    if [ -d /lustre/scratch123/tol/share/geval-data-prod/punchlist/${ref_name}_${dbVersion} ]
    then
        echo 'Directory Exists:\tcheck contents before -resume'
        exit 999
    else
        echo 'Creating Directory:\t/lustre/scratch123/tol/share/geval-data-prod/punchlist/${ref_name}_${dbVersion}'
        mkdir /lustre/scratch123/tol/share/geval-data-prod/punchlist/${ref_name}_${dbVersion}
        cp ${punchlist} /lustre/scratch123/tol/share/geval-data-prod/punchlist/${ref_name}_${dbVersion}/
    fi
    """
}
