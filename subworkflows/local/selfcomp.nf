#!/usr/bin/env nextflow

import java.math.RoundingMode;
import java.math.BigDecimal;

//
// MODULE IMPORT BLOCK
//
include { MUMMER                         } from '../../modules/nf-core/mummer/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UCSC_BEDTOBIGBED               } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                  } from '../../modules/nf-core/bedtools/sort/main'
include { SELFCOMP_SPLITFASTA            } from '../../modules/local/selfcomp_splitfasta'
include { SELFCOMP_MUMMER2BED            } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS                } from '../../modules/local/selfcomp_mapids'
include { CHUNKFASTA                     } from '../../modules/local/chunkfasta'
include { CAT_CAT                        } from '../../modules/nf-core/cat/cat/main'
include { SELFCOMP_ALIGNMENTBLOCKS       } from '../../modules/local/selfcomp_alignmentblocks'
include { CONCATBLOCKS                   } from '../../modules/local/concatblocks'
include { BEDTOOLS_MERGE                 } from '../../modules/nf-core/bedtools/merge/main'

workflow SELFCOMP {
    take:
    reference_tuple      // Channel: tuple [ val(meta), path(reference_file) ]
    dot_genome           // Channel: tuple [ val(meta), [ path(datafile) ] ]
    mummer_chunk         // Channel: val( int )
    motif_len            // Channel: val( int )
    selfcomp_as          // Channel: val( dot_as location )

    main:
    ch_versions             = Channel.empty()

    //
    // MODULE: SPLITS INPUT FASTA INTO 50KB WINDOWS
    //          EMITS A SINGLE FILE CONTAINING THESE WINDOWS
    //          THIS ACTS AS THE REFERENCE FOR GENOME.size() < 1GB
    //
    SELFCOMP_SPLITFASTA(
        reference_tuple
    )
    ch_versions             = ch_versions.mix( SELFCOMP_SPLITFASTA.out.versions )

    //
    // LOGIC: CALCULATE THE NUMBER OF GB WHICH WILL DICTATE THE NUMBER OF
    //          CHUNKS THE REFERENCE NEEDS TO BE SPLIT INTO
    //          ALSO CALCULATES THE NUMBER OF TOTAL WINDOWS NEEDED IN THE REFERENCE
    //
    reference_tuple
        .map{ it, file -> file.size()}
        .set { file_size }                  // Using set as TAP will force the pipeline to not complete successfully in some cases

    file_size
        .sum{it / 1e9}
        .collect { new BigDecimal (it).setScale(0, RoundingMode.UP) }
        .flatten()
        .set { chunk_number }

    //
    // MODULE: SPLIT REFERENCE FILE INTO 1GB CHUNKS
    //          THIS IS THE QUERY, AND REFERENCE IF GENOME.size() > 1GB
    //
    CHUNKFASTA(
        SELFCOMP_SPLITFASTA.out.fa,
        chunk_number
    )
    ch_versions         = ch_versions.mix( CHUNKFASTA.out.versions )

    //
    // LOGIC: STRIP META FROM QUERY, AND COMBINE WITH REFERENCE FILE
    //          THIS LEAVES US WITH n=( REFERENCE + QUERY) IF GENOME.SIZE() < 1GB
    //          OR n=((REFERENCE / 1E9) * (REFENCE / 1E9)) IF GENOME.SIZE() > 1GB
    //
    CHUNKFASTA.out.fasta
        .map{ meta, query ->
            query
        }
        .collect()                                              // Collect any output from CHUNKFASTA
        .map { it ->
            tuple(  [   len: it.size()   ],                     // Calc length of list
                    it
            )
        }
        .set { len_ch }                                         // tap out to preserve length of CHUNKFASTA list

    len_ch                                                      // tap swapped with set as tap stops pipeline completion
        .map { meta, files ->
            files
        }
        .flatten()                                              // flatten list into singles
        .combine(len_ch)                                        // re-add length information
        .combine(SELFCOMP_SPLITFASTA.out.fa)                    // add proposed reference, will be replaced by query list if > 1gb
        .map{                                                   // map all data together, if lenth of list was larger then 1
                                                                // indicating the original file was size() > 1Gb
            qry, len_meta, len_collected, ref_meta, ref ->
                tuple([ id: qry.toString().split('/')[-1],
                        sz: len_meta.len
                    ],
                    ( len_meta.len > 1 ? qry : ref ),            // Swap ref for query if list > 1
                    ( len_meta.len > 1 ? len_collected: [qry])   // Swap query for collected list of query if list > 1
                )
        }
        .transpose()                                             // Transpose the channel so that we have a channel for file in query
                                                                 // allows this to work on list of 1 and beyond
        .map { meta, ref, qry ->
            tuple(  [   id: meta.id,
                        sz: meta.sz,
                        it: qry.toString().split('/')[-1]        // get file name of the new query
                    ],
                    ref,
                    qry
            )
        }
        .set{ mummer_input }

    //
    // MODULE: ALIGNS 1GB CHUNKS TO 500KB CHUNKS
    //         EMITS MUMMER ALIGNMENT FILE
    //
    MUMMER(
        mummer_input
    )
    ch_versions             = ch_versions.mix( MUMMER.out.versions )

    //
    // LOGIC: COLLECT COORD FILES AND CONVERT TO LIST OF FILES
    //          ADD REFERENCE META
    //
    MUMMER.out.coords
        .map{ meta, file ->
            file
        }
        .collect()
        .toList()
        .combine( reference_tuple )
        .map { files, meta, ref ->
            tuple(  meta,
                    files
            )
        }
        .set { ch_mummer_files }

    //
    // MODULE: MERGES MUMMER ALIGNMENT FILES
    //
    CAT_CAT(
        ch_mummer_files
    )
    ch_versions             = ch_versions.mix( CAT_CAT.out.versions )

    //
    // MODULE: CONVERT THE MUMMER ALIGNMENTS INTO BED FORMAT
    //
    SELFCOMP_MUMMER2BED(
        CAT_CAT.out.file_out,
        motif_len
    )
    ch_versions             = ch_versions.mix( SELFCOMP_MUMMER2BED.out.versions )

    //
    // MODULE: GENERATE A LIST OF IDs AND GENOMIC POSITIONS OF SELFCOMPLEMENTARY REGIONS
    //         EMITS BED FILE
    //
    SELFCOMP_MAPIDS(
        SELFCOMP_MUMMER2BED.out.bedfile,
        SELFCOMP_SPLITFASTA.out.agp
    )
    ch_versions             = ch_versions.mix( SELFCOMP_MAPIDS.out.versions )

    //
    // MODULE: SORTS ABOVE OUTPUT BED FILE AND RETAINS BED SUFFIX
    //
    BEDTOOLS_SORT(
        SELFCOMP_MAPIDS.out.bedfile,
        []
    )
    ch_versions             = ch_versions.mix( BEDTOOLS_SORT.out.versions )

    //
    // MODULE: BUILD ALIGNMENT BLOCKS
    //
    SELFCOMP_ALIGNMENTBLOCKS(
        BEDTOOLS_SORT.out.sorted
    )
    ch_versions             = ch_versions.mix( SELFCOMP_ALIGNMENTBLOCKS.out.versions )

    //
    // MODULE: SORT BLOCKS FILES AND FILTER BY MOTIF LENGTH
    //
    CONCATBLOCKS(
        SELFCOMP_ALIGNMENTBLOCKS.out.blockfile
    )
    ch_versions             = ch_versions.mix( CONCATBLOCKS.out.versions )

    //
    // MODULE: CONVERTS ABOVE OUTPUT INTO BIGBED FORMAT
    //
    UCSC_BEDTOBIGBED(
        CONCATBLOCKS.out.chainfile,
        dot_genome.map{it[1]}, // Pulls file from tuple ( meta and file )
        selfcomp_as
    )
    ch_versions             = ch_versions.mix( UCSC_BEDTOBIGBED.out.versions )

    emit:
    ch_bigbed               = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions.ifEmpty(null)
}

