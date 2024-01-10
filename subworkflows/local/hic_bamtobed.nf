#!/usr/bin/env nextflow

// This subworkflow takes converts .bam to .bed for the hic_mapping subworkflow. 
// It runs markdup, sort and get paired contacts.
// Input - Assembled genomic fasta file, .bam file
// Output - sorted .bed and paired contact .bed

//
// MODULE IMPORT BLOCK
//
include { SAMTOOLS_MARKDUP                          } from '../../modules/nf-core/samtools/markdup/main'
include { BAMTOBED_SORT                             } from '../../modules/local/bamtobed_sort.nf'
include { GET_PAIRED_CONTACT_BED                    } from '../../modules/local/get_paired_contact_bed'


workflow HIC_BAMTOBED {
    take:
    bam_file            // Channel: tuple [ val(meta), path( file )      ]
    reference_tuple     // Channel: tuple [ val(meta), path( file )      ]

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: PREPARE MARKDUP INPUT
    //
    bam_file
        .combine( reference_tuple )
        .multiMap {  meta_bam, bam_file, meta_ref, ref ->
            bam      :   tuple(meta_bam, bam_file )
            reference      :   ref
        }
        .set { markdup_input }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MARKDUP (
        markdup_input.bam,
        markdup_input.reference
    )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MARKDUP.out.versions )

    //
    // MODULE: SAMTOOLS FILTER OUT DUPLICATE READS | BAMTOBED | SORT BED FILE
    //
    BAMTOBED_SORT(
        SAMTOOLS_MARKDUP.out.bam
    )
    ch_versions         = ch_versions.mix( BAMTOBED_SORT.out.versions )

    //
    // MODULE: GENERATE CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED( 
        BAMTOBED_SORT.out.sorted_bed 
    )
    ch_versions         = ch_versions.mix( GET_PAIRED_CONTACT_BED.out.versions )

    emit:
    paired_contacts_bed = GET_PAIRED_CONTACT_BED.out.bed
    sorted_bed          = BAMTOBED_SORT.out.sorted_bed
    versions            = ch_versions.ifEmpty(null)
}
