#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { MUMMER                         } from '../../modules/nf-core/mummer/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UCSC_BEDTOBIGBED               } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                  } from '../../modules/nf-core/bedtools/sort/main'
include { SELFCOMP_SPLITFASTA            } from '../../modules/local/selfcomp_splitfasta'
include { SELFCOMP_MUMMER2BED            } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS                } from '../../modules/local/selfcomp_mapids'
include { CHUNKFASTA                     } from '../../modules/local/chunkfasta'
include { CONCATMUMMER                   } from '../../modules/local/concatmummer'
include { SELFCOMP_ALIGNMENTBLOCKS       } from '../../modules/local/selfcomp_alignmentblocks'
include { SELFCOMP_MERGEBLOCKS           } from '../../modules/local/selfcomp_mergeblocks'
include { CONCATBLOCKS                   } from '../../modules/local/concatblocks'

workflow SELFCOMP {
    take:
        reference_tuple      // channel [id: sample_id], reference_file
        dot_genome           // Channel: [val(meta), [ datafile ]]
        mummer_chunk         // channel val( int )
        motif_len            // channel val( int )
        selfcomp_as          // channel val(dot_as location)

    main:
    ch_versions             = Channel.empty()
     
    // 
    // MODULE: SPLITS INPUT FASTA INTO 500KB CHUNKS
    //         EMITS CHUNKED FASTA
    //
    SELFCOMP_SPLITFASTA(reference_tuple)
    ch_versions             = ch_versions.mix(SELFCOMP_SPLITFASTA.out.versions)

    //
    // MODULE: SPLIT INPUT FASTA INTO 1GB CHUNKS
    //         EMITS CHUNKED FASTA
    //
    CHUNKFASTA(SELFCOMP_SPLITFASTA.out.fa, mummer_chunk)
    ch_versions             = ch_versions.mix(CHUNKFASTA.out.versions)

    //
    // LOGIC: CONVERTS ABOVE OUTPUTS INTO A SINGLE TUPLE
    //
    ch_query_tup = CHUNKFASTA.out.fas
        .map{
            meta, query -> 
            [query]
        }
        .flatten()

    ch_ref = SELFCOMP_SPLITFASTA.out.fa
        .map{
            meta, ref -> 
            ref
        }

    ch_mummer_input = ch_query_tup
        .combine(ch_ref)
        .map{ query, ref -> 
              tuple([id: query.toString().split('/')[-1] ], 
              ref, 
              query
              )
        }

    //
    // MODULE: ALIGNS 1GB CHUNKS TO 500KB CHUNKS
    //         EMITS MUMMER ALIGNMENT FILE
    //
    MUMMER( ch_mummer_input )
    ch_versions             = ch_versions.mix(MUMMER.out.versions)

    //
    // LOGIC: GROUPS OUTPUT INTO SINGLE TUPLE BASED ON REFERENCE META
    //
    MUMMER.out.coords
        .combine(reference_tuple)
        .map { coords_meta, coords, ref_meta, ref -> 
               tuple( ref_meta, 
               coords 
               ) 
        }
        .groupTuple(by:[0])
        .set{ ch_mummer_files }


    //
    // MODULE: MERGES MUMMER ALIGNMENT FILES
    //
    CONCATMUMMER(ch_mummer_files)
    ch_versions             = ch_versions.mix(CONCATMUMMER.out.versions)

    //
    // MODULE: CONVERT THE MUMMER ALIGNMENTS INTO BED FORMAT
    //
    SELFCOMP_MUMMER2BED(CONCATMUMMER.out.mummer, motif_len)
    ch_versions             = ch_versions.mix(SELFCOMP_MUMMER2BED.out.versions)

    //
    // MODULE: GENERATE A LIST OF IDs AND GENOMIC POSITIONS OF SELFCOMPLEMENTARY REGIONS
    //         EMITS BED FILE
    //
    SELFCOMP_MAPIDS(SELFCOMP_MUMMER2BED.out.bedfile, SELFCOMP_SPLITFASTA.out.agp)
    ch_versions             = ch_versions.mix(SELFCOMP_MAPIDS.out.versions)

    //
    // MODULE: SORTS ABOVE OUTPUT BED FILE AND RETAINS BED SUFFIX
    //
    BEDTOOLS_SORT(SELFCOMP_MAPIDS.out.bedfile, [])
    ch_versions             = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    //
    // MODULE: BUILD ALIGNMENT BLOCKS
    //
    SELFCOMP_ALIGNMENTBLOCKS(BEDTOOLS_SORT.out.sorted)
    ch_versions             = ch_versions.mix(SELFCOMP_ALIGNMENTBLOCKS.out.versions)

    //
    // LOGIC: CONVERTS ABOVE OUTPUTS INTO A LIST
    //
    SELFCOMP_ALIGNMENTBLOCKS.out.blockfile
        .map{ id, block -> 
              block
        }
        .flatten()
        .set{ch_blocks}
    
    //
    // LOGIC: CONVERTS CH_BLOCKS TO A LIST OF TUPLES
    //
    ch_blocks
        .map{ row ->
              tuple([id:row.toString().split('/')[-1]], 
              file(row)
            )
        }
        .set{ch_mergeblock_input}


    //
    // MODULE: BEDTOOLS MERGE ALIGNMENT BLOCKS
    //
    SELFCOMP_MERGEBLOCKS(ch_mergeblock_input)
    ch_versions             = ch_versions.mix(SELFCOMP_MERGEBLOCKS.out.versions)

    //
    // LOGIC: CONVERTS ABOVE OUTPUTS INTO A TUPLE
    //
    SELFCOMP_MERGEBLOCKS.out.mergedblocks
    .combine(reference_tuple)
    .map { merge_meta, mergedblocks, ref_meta, ref -> 
           tuple( ref_meta, 
           mergedblocks 
           ) 
         }
        .groupTuple(by:[0])
        .set{ch_merge}

    //
    // MODULE: MERGE ALL INDIVIDUAL BLOCKS FILES AND FILTER BY MOTIF LENGTH
    //
    CONCATBLOCKS(ch_merge)
    ch_versions             = ch_versions.mix(SELFCOMP_MERGEBLOCKS.out.versions)

    //
    // MODULE: CONVERTS ABOVE OUTPUT INTO BIGBED FORMAT
    //
    UCSC_BEDTOBIGBED(CONCATBLOCKS.out.chainfile, dot_genome.map{it[1]}, selfcomp_as)
    ch_versions             = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    ch_bigbed              = UCSC_BEDTOBIGBED.out.bigbed
    versions               = ch_versions.ifEmpty(null)
}
