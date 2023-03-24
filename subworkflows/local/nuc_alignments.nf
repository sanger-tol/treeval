include { MINIMAP2_ALIGN        } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE        } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/samtools/faidx/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/bedtools/sort/main'
include { BEDTOOLS_BAMTOBED     } from '../../modules/nf-core/bedtools/bamtobed/main'
include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/ucsc/bedtobigbed/main'


workflow NUC_ALIGNMENTS {
    take:
    reference_tuple
    reference_index
    nuc_files
    dot_genome
    intron_size

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: COLLECTION FROM GENE_ALIGNMENT IS A LIST OF ALL META AND ALL FILES
    //        BELOW CONVERTS INTO TUPLE FORMAT AND ADDS BOOLEANS FOR MINIMAP2_ALIGN
    //
    nuc_files
        .flatten()
        .buffer( size: 2 )
        .combine ( reference_tuple )
        .combine( intron_size )
        .map ( it ->
            tuple( [id:             it[0].id,
                    type:           it[0].type,
                    org:            it[0].org,
                    intron_size:    it[4],
                    single_end: true
                    ],
                    it[1],
                    it[3],
                    true,
                    false,
                    false
            )
        )
        .set { formatted_input }

    //
    // MODULE: ALIGNS REFERENCE FAIDX TO THE GENE_ALIGNMENT QUERY FILE FROM NUC_FILES
    //         EMITS ALIGNED BAM FILE
    //
    MINIMAP2_ALIGN (
        formatted_input.map { [it[0], it[1]] },
        formatted_input.map { it[2] },
        formatted_input.map { it[3] },
        formatted_input.map { it[4] },
        formatted_input.map { it[5] }
    )
    ch_versions     = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // LOGIC: CONVERTS THE MINIMAP OUTPUT TUPLE INTO A GROUPED TUPLE PER INPUT QUERY ORGANISM 
    //        AND DATA TYPE (RNA, CDS, DNA).
    //        EMITS THREE CHANNELS FOR THE GROUPED QUERY DATA REFERENCE AND REFERENCE INDEX
    //
    MINIMAP2_ALIGN.out.bam
        .map { meta, file ->
            tuple([id: meta.org, type: meta.type], file) } 
        .groupTuple( by: [0] )
        .combine( reference_tuple )
        .combine( reference_index )
        .multiMap { it ->
            nuc_grouped:    tuple( it[0], it[1] )
            reference:      it[-3]
            ref_index:      it[-1]
        }
        .set { merge_input }

    //
    // MODULE: MERGES THE BAM FILES FOUND IN THE GROUPED TUPLE IN REGARDS TO THE REFERENCE
    //         EMITS A MERGED BAM
    SAMTOOLS_MERGE (
        merge_input.nuc_grouped,
        merge_input.reference, 
        merge_input.ref_index
    )
    ch_versions     = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // MODULE: CONVERTS THE ABOVE MERGED BAM INTO BED FORMAT
    //
    BEDTOOLS_BAMTOBED { SAMTOOLS_MERGE.out.bam }

    //
    // MODULE: SORTS THE ABOVE BED FILE
    //
    BEDTOOLS_SORT ( BEDTOOLS_BAMTOBED.out.bed, [] )
    ch_versions     = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    //
    // LOGIC: COMBINES GENOME_FILE CHANNEL AND ABOVE OUTPUT, SPLITS INTO TWO CHANNELS
    //
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
    //
    // MODULE: CONVERTS GENOME FILE AND BED INTO A BIGBED FILE
    //
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
