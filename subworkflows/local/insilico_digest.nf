#!/usr/bin/env nextflow
//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

include { MAKECMAP_FA2CMAPMULTICOLOR } from '../../modules/sanger-tol/nf-core-modules/makecmap/fa2cmapmulticolor/main'
include { MAKECMAP_RENAMECMAPIDS } from '../../modules/sanger-tol/nf-core-modules/makecmap/renamecmapids/main'
nextflow.enable.dsl = 2

workflow INSILICO_DIGEST {

    main:

    input_fasta = [
        [ id: params.sample, single_end:false ], // meta map
        file(params.fasta, checkIfExists: true)
    ]

    MAKECMAP_FA2CMAPMULTICOLOR ( input_fasta, params.enzyme )
}