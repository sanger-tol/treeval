#!/usr/bin/env nextflow

//
// MODULE IMPORT BLOCK
//
include { CAT_CAT               } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK           } from "../../modules/nf-core/fastk/fastk/main"
include { FKUTILS_FKPROF        } from "../../modules/local/fkutils/fkprof/main"
include { GNU_SORT              } from "../../modules/nf-core/gnu/sort/main"
include { UCSC_BEDGRAPHTOBIGWIG } from "../../modules/nf-core/ucsc/bedgraphtobigwig/main"


workflow KMER_READ_COVERAGE {
    take:
    dot_genome                      // Channel: [ val(meta), path(file) ]
    reference                       // Channel: [ val(meta), path(file) ]
    read_ch                      // Channel: [ val(meta), path(file) ]
    kmer_prof_file                  // Channel: [ val( meta[id, kmer] ), path(file) ]

    main:
    ch_versions             = Channel.empty()

    kmer_len    = kmer_prof_file.map( it -> it[0].kmer ).toInteger()
    kmer_file   = kmer_prof_file.map( it -> it[1] )

    if ( kmer_file.toString().endsWith('.ktab')) {

        //
        // LOGIC: GET FOLDER OF THE FASTK DATA - CONTENTS ARE REQUIRED FOR THE RUNNING OF THE FKPROF
        //          EVEN THOUGH IT ONLY WANTS THE KTAB FILE PASSED IN
        //
        kmer_prof_file
            .map { meta, file ->
                tuple ( meta,
                        file.getParent().toString())
            }
        .set { corrected_input }

        //
        // MODULE: PROFILE THE KMER SPECTRA
        //
        FKUTILS_FKPROF(
            reference,
            corrected_input
        )
        ch_versions             = ch_versions.mix( FKUTILS_FKPROF.out.versions )

    } else {

        //
        // MODULE: GETS PACBIO READ PATHS FROM READS_PATH
        //
        ch_grabbed_read_paths   = GrabFiles( read_ch )

        //
        // MODULE: JOIN PACBIO READ
        //
        CAT_CAT(
            ch_grabbed_read_paths
        )
        ch_versions             = ch_versions.mix( CAT_CAT.out.versions )

        //
        // LOGIC: ADDING THE KMER LENGTH INTO THE CHANNEL
        //
        CAT_CAT.out.file_out
            .combine( kmer_len )
            .map { meta, file, kmer_size ->
                tuple ( [   id  : meta.id,
                            kmer: kmer_size ],
                        file
                )
            }
            .set { cat_and_kmer }

        //
        // MODULE: COUNT KMERS BASED ON KMER_PROFILE meta.kmer
        //
        FASTK_FASTK(
            cat_and_kmer
        )
        ch_versions             = ch_versions.mix( FASTK_FASTK.out.versions )

        //
        // LOGIC: GET FOLDER OF THE FASTK DATA - CONTENTS ARE REQUIRED FOR THE RUNNING OF THE FKPROF
        //          EVEN THOUGH IT ONLY WANTS THE KTAB FILE PASSED IN
        //
        FASTK_FASTK.out.ktab
            .map { meta, files ->
                tuple(  meta,
                        files[0].getParent().toString()
                )
            }
        .set { filtered_ktabs }

        //
        // MODULE: PROFILE THE KMER SPECTRA
        //
        FKUTILS_FKPROF(
            reference,
            filtered_ktabs
        )
        ch_versions             = ch_versions.mix( FKUTILS_FKPROF.out.versions )
    }

    //
    // MODULE: SORT THE FKUTILS/FKPROF BED FILE
    //
    GNU_SORT (
        FKUTILS_FKPROF.out.bed
    )
    ch_versions                 = ch_versions.mix( GNU_SORT.out.versions )

    //
    // MODULE: CONVERT BED TO BIGWIG
    // LOGIC: STRIP THE META FROM THE GENOME FILE FOR THE SORTING TO WORK
    //
    genome_file = dot_genome.map { it -> it[1] }

    GNU_SORT.out.sorted
        .combine( kmer_len )
        .map { meta, file, kmer ->
            tuple ( [   id:     meta.id,
                        kmer:   kmer    ],
                    file)
        }
        .set {kmer_and_file}

    UCSC_BEDGRAPHTOBIGWIG(
        kmer_and_file,
        genome_file
    )
    ch_versions                 = ch_versions.mix( UCSC_BEDGRAPHTOBIGWIG.out.versions )


    emit:
    bigwig                      = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // meta, bigwig
    versions                    = ch_versions.ifEmpty( null )

}

process GrabFiles {
    label 'process_tiny'

    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.{fa,fasta}.{gz}")

    "true"
}
