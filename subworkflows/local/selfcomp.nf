#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { MUMMER                } from '../../modules/nf-core/modules/mummer/main'
include { SAMTOOLS_FAIDX        } from '../../modules/nf-core/modules/samtools/faidx/main'
include { SELFCOMP_SPLITFASTA   } from '../../modules/sanger-tol/nf-core-modules/selfcomp/splitfasta/main'
include { SELFCOMP_MUMMER2BED   } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS       } from '../../modules/local/selfcomp_mapids'
include { CHUNKFASTA            } from '../../modules/local/chunkfasta'
include { CONCATMUMMER          } from '../../modules/local/concatmummer'
include { UCSC_BEDTOBIGBED      } from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT         } from '../../modules/nf-core/modules/bedtools/sort/main'

workflow SELFCOMP {
    take:
        reference_tuple     // channel [id: sample_id], reference_file
        dot_genome          // Channel: [val(meta), [ datafile ]]
        mummer_chunk        // channel val( int )
        motif_len           // channel val( int )
        selfcomp_as         // channel val(dot_as location)

    main:
    ch_versions = Channel.empty()
     
    // Split fasta sequences into 500kb length fasta sequences
    SELFCOMP_SPLITFASTA(reference_tuple)
    ch_versions = ch_versions.mix(SELFCOMP_SPLITFASTA.out.versions)

    // Chunk the fasta file for processing (1Gb chunks)
    CHUNKFASTA(SELFCOMP_SPLITFASTA.out.fa, mummer_chunk)
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
    MUMMER.out.coords
        .combine(reference_tuple)
        .map { coords_meta, coords, ref_meta, ref -> tuple(ref_meta, coords) }
        .groupTuple(by:[0])
        .set{ ch_mummer_files }

    CONCATMUMMER(ch_mummer_files)
    ch_versions = ch_versions.mix(CONCATMUMMER.out.versions)

    SELFCOMP_MUMMER2BED(CONCATMUMMER.out.mummer, motif_len)
    ch_versions = ch_versions.mix(SELFCOMP_MUMMER2BED.out.versions)

    SELFCOMP_MAPIDS(SELFCOMP_MUMMER2BED.out.bedfile, SELFCOMP_SPLITFASTA.out.agp)
    ch_versions = ch_versions.mix(SELFCOMP_MAPIDS.out.versions)

    BEDTOOLS_SORT(SELFCOMP_MAPIDS.out.bedfile, "bed")
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    UCSC_BEDTOBIGBED(BEDTOOLS_SORT.out.sorted, dot_genome.map{it[1]}, selfcomp_as)
    ch_bigbed = UCSC_BEDTOBIGBED.out.bigbed
    ch_versions = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_bigbed
    versions = ch_versions
}