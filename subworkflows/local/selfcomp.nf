#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { MUMMER } from '../../modules/nf-core/modules/mummer/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'
include { SELFCOMP_SPLITFASTA as SPLITFASTA } from '../../modules/sanger-tol/nf-core-modules/selfcomp/splitfasta/main'
include { SELFCOMP_MUMMER2BED as MUMMER2BED } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mummer2bed/main'
include { SELFCOMP_MAPIDS as MAPIDS } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mapids/main'
include { CHUNKFASTA } from '../../modules/local/chunkfasta'
include { CONCATMUMMER } from '../../modules/local/concatmummer'
include { UCSC_BEDTOBIGBED} from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'

workflow SELFCOMP {

    main:
    ch_versions = Channel.empty()

    // Inputs
    input_fasta = params.reference
    sample_name = params.assembly.sample
    genome_size = params.genome_size
    number_of_chunks = 6
    motiflen = 0
    
    // Split fasta sequences into 500kb length fasta sequences
    SPLITFASTA([[id:sample_name, single_end: true], input_fasta])
    ch_split_fa = SPLITFASTA.out.fa
    ch_split_agp = SPLITFASTA.out.agp

    // Chunk the fasta file for processing (1Gb chunks)
    CHUNKFASTA(ch_split_fa, number_of_chunks)
    ch_chunked_fas = CHUNKFASTA.out.fas

    // Run mummer on each chunk, generating .coords files for each.
    ch_query_tup = ch_chunked_fas
    .map{
        meta, query -> [query]
        }
    .flatten()
    ch_mummer_input = ch_query_tup
    .map{
        query -> tuple([id: query.toString().split('/')[-1] ], file(input_fasta), query)
    }
    MUMMER( ch_mummer_input )
    ch_mummer_out = MUMMER.out.coords

    // Concatenate mummer files.
    ch_mummer_files = ch_mummer_out
        .map{
        meta, query -> [query]
        }
        .collect()

    ch_concatmummer_input = ch_mummer_files
        .map{
            query -> tuple([id: sample_name], query)
        }
        .view()

    CONCATMUMMER( ch_concatmummer_input )
    ch_concat_mummer_out = CONCATMUMMER.out.mummer

    // Run mummer2bed module (.mummer, motiflen) -> bed.
    MUMMER2BED(ch_concat_mummer_out, motiflen)
    ch_bed = MUMMER2BED.out.bedfile

    // Run mapids module (.bed, .agp) -> .bed
    MAPIDS(ch_bed, ch_split_agp)
    ch_mapped_bed = MAPIDS.out.bedfile

    // Convert bed to bigbed, requires chromsize and assembly.as
    UCSC_BEDTOBIGBED(ch_mapped_bed, genome_size)
    ch_bigbed = UCSC_BEDTOBIGBED.out.bigbed

    emit:
    ch_versions

}