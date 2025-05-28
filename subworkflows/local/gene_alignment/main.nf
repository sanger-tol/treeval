#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

//
// SUBWORKFLOW IMPORT BLOCK
//
include { PEP_ALIGNMENTS                    } from '../pep_alignments/main'
include { NUC_ALIGNMENTS as GEN_ALIGNMENTS  } from '../nuc_alignments/main'
include { NUC_ALIGNMENTS as RNA_ALIGNMENTS  } from '../nuc_alignments/main'
include { NUC_ALIGNMENTS as CDS_ALIGNMENTS  } from '../nuc_alignments/main'

workflow GENE_ALIGNMENT {
    take:
    dot_genome          // Channel: [ val(meta), path(file) ]
    reference_tuple     // Channel: [ val(meta), path(file) ]
    reference_index     // Channel: [ val(meta), path(file) ]
    alignment_genesets  // Channel: val(geneset_id)
    intron_size         // Channel: val(50k)
    as_files            // Channel: [ val(meta), path(file) ]

    main:
    ch_versions         = Channel.empty()

    reference_tuple
        .map{ meta, file ->
            "${meta.class}"
        }
        .set { assembly_class }


    //
    // LOGIC: TAKES A SINGLE LIKE CSV STRING AND CONVERTS TO LIST OF VALUES
    //          LIST IS MERGED WITH DATA_DIRECTORY AND ORGANISM_CLASS
    //
    ch_data             = alignment_genesets
                            .splitCsv()
                            .flatten()

    //
    // LOGIC:   COMBINE CH_DATA WITH ALIGNMENT_DIR AND ASSEMBLY_CLASS
    //          CONVERTS THESE VALUES INTO A PATH AND DOWNLOADS IT, THEN TURNS IT TO A TUPLE OF
    //          [ [ META.ID, META.TYPE, META.ORG ], GENE_ALIGNMENT_FILE ]
    //          DATA IS THEN BRANCHED BASED ON META.TYPE TO THE APPROPRIATE
    //          SUBWORKFLOW
    //
    ch_data
        .map {
            geneset_path ->
                file(geneset_path)
        }
        .splitCsv( header: true, sep:',')
        .map{ row ->
            tuple([ org:    row.org,
                    type:   row.type,
                    id:     row.data_file.split('/')[-1].split('.MOD.')[0]
                ],
                file(row.data_file)
            )}
        .branch {
            pep: it[0].type  == 'pep'
            gen: it[0].type  == 'cdna'
            rna: it[0].type  == 'rna'
            cds: it[0].type  == 'cds'
        }
        .set {ch_alignment_data}

    pep_files = ch_alignment_data.pep.collect()
    gen_files = ch_alignment_data.gen.collect()
    rna_files = ch_alignment_data.rna.collect()
    cds_files = ch_alignment_data.cds.collect()

    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR PEPTIDE DATA, EMITS GFF AND TBI
    //
    PEP_ALIGNMENTS (    reference_tuple,
                        pep_files,
    )
    ch_versions = ch_versions.mix(PEP_ALIGNMENTS.out.versions)

    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR RNA, NUCLEAR AND COMPLEMENT_DNA DATA, EMITS BIGBED
    //
    GEN_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        gen_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( GEN_ALIGNMENTS.out.versions )

    CDS_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        cds_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( CDS_ALIGNMENTS.out.versions )

    RNA_ALIGNMENTS (    reference_tuple,
                        reference_index,
                        rna_files,
                        dot_genome,
                        intron_size
    )
    ch_versions = ch_versions.mix( RNA_ALIGNMENTS.out.versions )

    emit:
    pep_gff             = PEP_ALIGNMENTS.out.tbi_gff
    gff_file            = PEP_ALIGNMENTS.out.gff_file
    gen_bb_files        = GEN_ALIGNMENTS.out.nuc_alignment
    rna_bb_files        = RNA_ALIGNMENTS.out.nuc_alignment
    cds_bb_files        = CDS_ALIGNMENTS.out.nuc_alignment
    versions            = ch_versions
}
