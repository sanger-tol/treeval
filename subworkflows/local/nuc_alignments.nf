include { MINIMAP2_ALIGN        } from '../../modules/nf-core/modules/minimap2/align/main.nf'
include { SAMTOOLS_MERGE        } from '../../modules/nf-core/modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/modules/bedtools/sort/main'
include { BEDTOOLS_BAMTOBED     } from '../../modules/sanger-tol/nf-core-modules/bedtools/bamtobed/main'
include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'


workflow NUC_ALIGNMENTS {
    take:
    reference_tuple
    nuc_files
    dot_genome
    intron_size

    main:
    ch_versions         = Channel.empty()

    nuc_files
        .flatten()
        .buffer( size: 2 )
        .combine ( reference_tuple )
        .combine( intron_size )
        .map ( it ->
            tuple( [id:         it[0].id,
                    type:       it[0].type,
                    org:        it[0].org,
                    single_end: true
                    ],
                    it[1],
                    it[3],
                    true,
                    false,
                    false,
                    it[4]
            )
        )
        .set { formatted_input }

    SAMTOOLS_FAIDX ( reference_tuple )
    ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    MINIMAP2_ALIGN (
        formatted_input.map { [it[0], it[1]] },
        formatted_input.map { it[2] },
        formatted_input.map { it[3] },
        formatted_input.map { it[4] },
        formatted_input.map { it[5] },
        formatted_input.map { it[6] }
    )
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    MINIMAP2_ALIGN.out.bam
        .map { meta, file ->
            tuple([id: meta.org, type: meta.type], file) } 
        .groupTuple( by: [0] )
        .combine( reference_tuple )
        .combine( SAMTOOLS_FAIDX.out.fai )
        .multiMap { it ->
            nuc_grouped:    tuple( it[0], it[1] )
            reference:      it[-3]
            ref_index:      it[-1]
        }
        .set { merge_input }

    SAMTOOLS_MERGE (
        merge_input.nuc_grouped,
        merge_input.reference, 
        merge_input.ref_index
    )
    ch_versions     = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    BEDTOOLS_BAMTOBED { SAMTOOLS_MERGE.out.bam }

    BEDTOOLS_SORT ( BEDTOOLS_BAMTOBED.out.bed, 'sorted.bed' )
    ch_versions     = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    BEDTOOLS_SORT.out.sorted
        .combine( dot_genome )
        .multiMap { it ->
            bed_file:   tuple( [    id:     it[0].id,
                                    type:   it[0].type
                                ],
                                it[1] )
            dot_genome: it[3]    
        }
        .set { ucsc_input }

    ucsc_input.bed_file.view()

    UCSC_BEDTOBIGBED (
        ucsc_input.bed_file,
        ucsc_input.dot_genome,
        []
    )

    ch_versions     = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )

    emit:
    nuc_alignment   = UCSC_BEDTOBIGBED.out.bigbed.collect()
    versions        = ch_versions.ifEmpty(null)
}