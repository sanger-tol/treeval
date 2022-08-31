#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { MUMMER } from '../../modules/nf-core/modules/mummer/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'
include { SELFCOMP_SPLITFASTA } from '../../modules/sanger-tol/nf-core-modules/selfcomp/splitfasta/main'
include { SELFCOMP_MUMMER2BED } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mummer2bed/main'
include { SELFCOMP_MAPIDS } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mapids/main'
include { CHUNKFASTA } from '../../modules/local/chunkfasta'
include { CONCATMUMMER } from '../../modules/local/concatmummer'
include { UCSC_BEDTOBIGBED} from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT } from '../../modules/nf-core/modules/bedtools/sort/main'

workflow SELFCOMP {

    main:
    ch_versions = Channel.empty()

    // Inputs
    input_fasta = params.reference
    sample_name = params.assembly.sample
    genome_size = params.genome_size
    number_of_chunks = params.self_comp.mummer_chunk
    motiflen = params.self_comp.motif_len
    
    // Split fasta sequences into 500kb length fasta sequences
    SELFCOMP_SPLITFASTA([[id:sample_name, single_end: true], input_fasta])
    ch_versions = ch_versions.mix(SELFCOMP_SPLITFASTA.out.versions)

    // Chunk the fasta file for processing (1Gb chunks)
    CHUNKFASTA(SELFCOMP_SPLITFASTA.out.fa, number_of_chunks)
    ch_versions = ch_versions.mix(CHUNKFASTA.out.versions)

    // Run mummer on each chunk, generating .coords files for each.
    ch_query_tup = CHUNKFASTA.out.fas
        .map{
            meta, query -> [query]
        }
        .flatten()

    ch_ref = SELFCOMP_SPLITFASTA.out.fa
        .map{
            meta, ref -> ref
        }

    ch_mummer_input = ch_query_tup
        .combine(ch_ref)
        .map{
            query, ref -> tuple([id: query.toString().split('/')[-1] ], ref, query)
        }

    MUMMER( ch_mummer_input )
    ch_versions = ch_versions.mix(MUMMER.out.versions)

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
    ch_versions = ch_versions.mix(CONCATMUMMER.out.versions)

    SELFCOMP_MUMMER2BED(CONCATMUMMER.out.mummer, motiflen)
    ch_versions = ch_versions.mix(SELFCOMP_MUMMER2BED.out.versions)

    SELFCOMP_MAPIDS(SELFCOMP_MUMMER2BED.out.bedfile, SPLITFASTA.out.agp)
    ch_versions = ch_versions.mix(SELFCOMP_MAPIDS.out.versions)

    BEDTOOLS_SORT(SELFCOMP_MAPIDS.out.bedfile, "bed")
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    UCSC_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, genome_size)
    ch_bigbed = UCSC_BEDTOBIGBED.out.bigbed
    ch_versions = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_bigbed
    ch_versions
}