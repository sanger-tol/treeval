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
include { UCSC_BEDTOBIGBED as SELFCOMP_BEDTOBIGBED} from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/modules/bedtools/sort/main'

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

    // Chunk the fasta file for processing (1Gb chunks)
    CHUNKFASTA(SPLITFASTA.out.fa, number_of_chunks)

    // Run mummer on each chunk, generating .coords files for each.
    ch_query_tup = CHUNKFASTA.out.fas
        .map{
            meta, query -> [query]
        }
        .flatten()

    ch_ref = SPLITFASTA.out.fa
        .map{
            meta, ref -> ref
        }

    ch_mummer_input = ch_query_tup
        .combine(ch_ref)
        .map{
            query, ref -> tuple([id: query.toString().split('/')[-1] ], ref, query)
        }

    MUMMER( ch_mummer_input )

    // Concatenate mummer files.
    ch_mummer_files = MUMMER.out.coords
        .map{
            meta, query -> query
        }
        .collect()
        .map{
            query -> tuple([id: sample_name], query)
        }

    CONCATMUMMER(ch_mummer_files)

    MUMMER2BED(CONCATMUMMER.out.mummer, motiflen)

    MAPIDS(MUMMER2BED.out.bedfile, SPLITFASTA.out.agp)

    BEDTOOLS_SORT(MAPIDS.out.bedfile, "bed")

    SELFCOMP_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, genome_size)
    ch_bigbed = SELFCOMP_BEDTOBIGBED.out.bigbed

    emit:
    ch_versions
}