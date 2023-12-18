include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SE        } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE                             } from '../../modules/nf-core/samtools/merge/main'

workflow SE_MAPPING {

    take:
    reference_tuple          // Channel [ val(meta), path(file) ]
    assembly_path            // Channel path(file)
    pacbio_tuple             // Channel [ val(meta), path(file) ]
    reads_type                // Channel val( str )

    main:
    ch_versions     = Channel.empty()
    ch_align_bams   = Channel.empty()

    //
    // PROCESS: GETS PACBIO READ PATHS FROM READS_PATH
    //
    ch_grabbed_reads_path       = GrabFiles( pacbio_tuple )

    ch_grabbed_reads_path
        .map { meta, files ->
            tuple( files )
        }
        .flatten()
        .set { ch_reads_path }

    //
    // PROCESS: MAKE MINIMAP INPUT CHANNEL AND MAKE BRANCHES BASED ON INPUT READ TYPE
    //
    reference_tuple
        .combine( ch_reads_path )
        .combine( reads_type )
        .map { meta, ref, reads_path, reads_type ->
            tuple(
                [   id          : meta.id,
                    single_end  : true,
                    readtype    : reads_type.toString()
                ],
                reads_path,
                ref,
                true,
                false,
                false,
                reads_type
            )
        }
        .set { minimap_se_input }

    //
    // PROCESS: MULTIMAP TO MAKE BOOLEAN ARGUMENTS FOR MINIMAP HIFI MAPPING INPUT
    //
    minimap_se_input
        .multiMap { meta, reads_path, ref, bam_output, cigar_paf, cigar_bam, reads_type ->
            read_tuple          : tuple( meta, reads_path)
            ref                 : ref
            bool_bam_ouput      : bam_output
            bool_cigar_paf      : cigar_paf
            bool_cigar_bam      : cigar_bam
        }
        .set { se_input }

    //
    // MOUDLES: MAPPING DIFFERENT TYPE OF READ AGAINIST REFERENCE
    //

    MINIMAP2_ALIGN_SE (
            se_input.read_tuple,
            se_input.ref,
            se_input.bool_bam_ouput,
            se_input.bool_cigar_paf,
            se_input.bool_cigar_bam
    )
    ch_bams = MINIMAP2_ALIGN_SE.out.bam
   

    ch_bams
        .map { meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [ id    : file[0].toString().split('/')[-1].split('_')[0] ], // Change sample ID
                file
            )
        }
        .set { collected_files_for_merge }

    //
    // MODULE: MERGE ALL OUTPUT BAM
    //
    SAMTOOLS_MERGE(
        collected_files_for_merge,
        reference_tuple,
        [[],[]]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    versions       = ch_versions.ifEmpty(null)
    mapped_bam     = SAMTOOLS_MERGE.out.bam
}

process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.{fa,fasta}.{gz}")

    "true"
}
