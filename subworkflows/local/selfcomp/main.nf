#!/usr/bin/env nextflow



//
// MODULE IMPORT BLOCK
//
include { MUMMER                                 } from '../../modules/nf-core/mummer/main'
include { UCSC_BEDTOBIGBED                       } from '../../modules/nf-core/ucsc/bedtobigbed/main'
include { BEDTOOLS_SORT                          } from '../../modules/nf-core/bedtools/sort/main'
include { SELFCOMP_SPLITFASTA                    } from '../../modules/local/selfcomp_splitfasta'
include { SELFCOMP_MUMMER2BED                    } from '../../modules/local/selfcomp_mummer2bed'
include { SELFCOMP_MAPIDS                        } from '../../modules/local/selfcomp_mapids'
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_REF      } from '../../modules/nf-core/seqkit/split2/main'
include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_QUERY    } from '../../modules/nf-core/seqkit/split2/main'
include { CAT_CAT                                } from '../../modules/nf-core/cat/cat/main'
include { SELFCOMP_ALIGNMENTBLOCKS               } from '../../modules/local/selfcomp/alignmentblocks/main'
include { CONCAT_BLOCKS                           } from '../../modules/local/concat/blocks/main'

/*def processPaths(mytuple, prefix) {
    if (mytuple == null || mytuple.isEmpty() || mytuple[0] == null) {
        println "ERROR: processPaths received an empty or null tuple"
        return []
    }

    def pathList = mytuple[0][1] ?: []  // Safeguard against null
    def result = []

    pathList.eachWithIndex { pathString, idx ->
        def idxStr = "${prefix}${idx + 1}"
        result.add([[id: idxStr], pathString])
    }

    return result
}*/


workflow SELFCOMP {
    take:
    reference_tuple      // Channel: tuple [ val(meta), path(reference_file) ]
    dot_genome           // Channel: tuple [ val(meta), [ path(datafile) ] ]
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

    SELFCOMP_SPLITFASTA.out.fa
    .combine ( ref_chunk_number )
    .map{ meta, fastaFile, chunk_number ->
        tuple( [id:             meta.id,
                file_type:      "fasta",
                single_end:     true,
                cn:    chunk_number
                ],
                fastaFile
            )
    }
    .set { windowed_fasta_ref_ch }


    //
    // LOGIC: QUERY CHUNKS SHOULD BE UNDER 0.5GB PER CHUNK
    //
    reference_tuple
    .map { it, file ->
            def sizeInGB = (file.size() / 1_073_741_824.0)  / 0.5
            sizeInGB < 1 ? 1 : sizeInGB.toInteger()  // Conditional operator for the logic
    }
    .set { query_chunk_number }

    SELFCOMP_SPLITFASTA.out.fa
    .combine ( query_chunk_number )
    .map{ meta, fastaFile, chunk_number ->
        tuple( [id:             meta.id,
                file_type:      "fasta",
                single_end:     true,
                cn:    chunk_number
                ],
                fastaFile
            )
    }
    .set { windowed_fasta_query_ch }


    //
    // MODULE: SPLIT QUERY FILE INTO 1GB CHUNKS
    //          THIS IS THE QUERY, AND REFERENCE IF GENOME.size() > 1GB
    //
    SEQKIT_SPLIT_QUERY(
        windowed_fasta_query_ch
    )
    ch_versions         = ch_versions.mix(SEQKIT_SPLIT_QUERY.out.versions)

    SEQKIT_SPLIT_REF(
        windowed_fasta_ref_ch
    )
    ch_versions         = ch_versions.mix(SEQKIT_SPLIT_REF.out.versions)

    SEQKIT_SPLIT_REF.out.reads
    .map { meta, myfiles ->
        myfiles
    }
    .flatMap { it -> it}
    .set { ref_chunks }

    SEQKIT_SPLIT_QUERY.out.reads
    .map { meta, myfiles ->
        myfiles
    }
    .flatMap { it -> it}
    .set { query_chunks }

    //
    // LOGIC: CONSTRUCT MUMMER INPUT CHANNEL
    //
    ref_chunks
    .combine(query_chunks)
    .map { ref, q ->
        tuple([id: "${file(ref).getBaseName()}_${file(q).getBaseName()}"], file(ref), file(q))
    }
    .set { mummer_input }


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
        CAT_CAT.out.file_out
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
    ch_bigbed               = UCSC_BEDTOBIGBED.out.bigbed
    versions                = ch_versions.ifEmpty(null)
    }
