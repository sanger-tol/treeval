#!/usr/bin/env nextflow

import java.math.RoundingMode;
import java.math.BigDecimal;

//
// MODULE IMPORT BLOCK
//
include { MUMMER                                 } from '../../modules/nf-core/mummer/main'
include { UCSC_BEDTOBIGBED                       } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                          } from '../../modules/nf-core/bedtools/sort/main'
include { SELFCOMP_SPLITFASTA                    } from '../../modules/local/selfcomp_splitfasta'
include { SELFCOMP_MUMMER2BED                    } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS                        } from '../../modules/local/selfcomp_mapids'
include { SEQKIT_SPLIT as SEQKIT_SPLIT_REF       } from '../../modules/local/seqkit/split/main'
include { SEQKIT_SPLIT as SEQKIT_SPLIT_QUERY     } from '../../modules/local/seqkit/split/main'
include { CAT_CAT                                } from '../../modules/nf-core/cat/cat/main'
include { SELFCOMP_ALIGNMENTBLOCKS               } from '../../modules/local/selfcomp_alignmentblocks'
include { CONCATBLOCKS                           } from '../../modules/local/concatblocks'

workflow SELFCOMP {
    take:
    reference_tuple      // Channel: tuple [ val(meta), path(reference_file) ]
    dot_genome           // Channel: tuple [ val(meta), [ path(datafile) ] ]
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
    // LOGIC: REFERENCE SHOULD BE UNDER 1GB TO OPTIMIZE MEMORY USAGE
    //
    reference_tuple
    .map { it, file ->
           def sizeInGB = (file.size() / 1_073_741_824.0) + 0.5
           sizeInGB < 1 ? 1 : sizeInGB.toInteger()  // Conditional operator for the logic
    }
    .set { ref_chunk_number }
    
    //
    // LOGIC: QUERY CHUNKS SHOULD BE UNDER 0.5GB PER CHUNK
    //
    reference_tuple
    .map { it, file ->
           def sizeInGB = (file.size() / 1_073_741_824.0)  / 0.5
           sizeInGB < 1 ? 1 : sizeInGB.toInteger()  // Conditional operator for the logic
    }
    .set { query_chunk_number }


    //
    // MODULE: SPLIT QUERY FILE INTO 1GB CHUNKS
    //          THIS IS THE QUERY, AND REFERENCE IF GENOME.size() > 1GB
    //
    SEQKIT_SPLIT_QUERY(
        SELFCOMP_SPLITFASTA.out.fa,
        query_chunk_number
    )
    ch_versions         = ch_versions.mix(SEQKIT_SPLIT_QUERY.out.versions)

    SEQKIT_SPLIT_REF(
        SELFCOMP_SPLITFASTA.out.fa,
        ref_chunk_number
    )
    ch_versions         = ch_versions.mix(SEQKIT_SPLIT_REF.out.versions)

    //
    // LOGIC: RECONSTRUCT QUERY TUPLE
    //
    SEQKIT_SPLIT_QUERY.out.fasta
    .toList()
    .collect { tuple ->
  
            def metadata = tuple[0][0]  
            def paths = tuple[0][1]
            println paths     
            def pathList = paths.collect()
                
            def result = []
            pathList.eachWithIndex { pathString, idx ->
                    def qIdx = "query_${idx + 1}"
                    result.add([[id:qIdx], pathString])
                }
                
            return result
          
    }
    .flatMap { it -> it }
    .set { query_chunks }

    //
    // LOGIC: RECONSTRUCT REFERENCE TUPLE
    //

    SEQKIT_SPLIT_REF.out.fasta
    .toList()
    .collect { tuple ->
  
            def metadata = tuple[0][0]  
            def paths = tuple[0][1]
            println paths     
            def pathList = paths.collect()
                
            def result = []
            pathList.eachWithIndex { pathString, idx ->
                    def rIdx = "ref_${idx + 1}"
                    result.add([[id:rIdx], pathString])
                }
                
            return result
          
    }
    .flatMap { it -> it }
    .set { ref_chunks }

    //ref_chunks.view()
    //query_chunks.view()
    //
    // LOGIC: CONSTRUCT MUMMER INPUT CHANNEL
    //
    ref_chunks
        .combine(query_chunks)
        .map { refID, refpath, queryID, qpath -> 
            def myid = "${refID.id}_${queryID.id}"
            tuple ( [ id  : myid,
                      rid : refID.id,
                      qid : queryID.id
                    ], 
                      refpath, 
                      qpath 
                  )
        } 
        .set { mummer_input }
    
    mummer_input.view()
    
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
    // LOGIC: ADDING LINE COUNT TO THE FILE FOR BETTER RESOURCE USAGE
    //
    SELFCOMP_MAPIDS.out.bedfile
        .map { meta, file ->
            tuple ( [   id:     meta.id,
                        lines:  file.countLines()
                    ],
                    file
            )
        }
        .set { bedtools_input }

    //
    // MODULE: SORTS ABOVE OUTPUT BED FILE AND RETAINS BED SUFFIX
    //
    BEDTOOLS_SORT(
        bedtools_input,
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
    //ch_bigbed               = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions.ifEmpty(null)
    }
