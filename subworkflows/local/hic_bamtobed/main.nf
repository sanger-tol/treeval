#!/usr/bin/env nextflow

// This subworkflow takes converts .bam to .bed for the hic_mapping subworkflow.
// It runs markdup, sort and get paired contacts.
// Input - Assembled genomic fasta file, .bam file
// Output - sorted .bed and paired contact .bed

//
// MODULE IMPORT BLOCK
//
include { BEDTOOLS_BAMTOBEDSORT  } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { GET_PAIRED_CONTACT_BED } from '../../../modules/local/get/paired_contact_bed/main'


workflow HIC_BAMTOBED {
    take:
    ch_bam_file // Channel: tuple [ val(meta), path( file )      ]

    main:
    ch_versions = channel.empty()

    //
    // MODULE: SAMTOOLS FILTER OUT DUPLICATE READS | BAMTOBED | SORT BED FILE
    //
    BEDTOOLS_BAMTOBEDSORT(ch_bam_file)

    //
    // MODULE: GENERATE CONTACT PAIRS
    //
    GET_PAIRED_CONTACT_BED(
        BEDTOOLS_BAMTOBEDSORT.out.sorted_bed
    )
    ch_versions = ch_versions.mix(GET_PAIRED_CONTACT_BED.out.versions)

    emit:
    paired_contacts_bed = GET_PAIRED_CONTACT_BED.out.bed
    sorted_bed          = BEDTOOLS_BAMTOBEDSORT.out.sorted_bed
    versions            = ch_versions
}
