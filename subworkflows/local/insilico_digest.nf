#!/usr/bin/env nextflow
//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

include { MAKECMAP_FA2CMAPMULTICOLOR } from '../../modules/sanger-tol/nf-core-modules/makecmap/fa2cmapmulticolor/main'
include { MAKECMAP_RENAMECMAPIDS } from '../../modules/sanger-tol/nf-core-modules/makecmap/renamecmapids/main'
include { MAKECMAP_CMAP2BED } from '../modules/sanger-tol/nf-core-modules/makecmap/cmap2bed/main’


nextflow.enable.dsl = 2

workflow INSILICO_DIGEST {

    main:

    myid = params.sample + "_" + params.enzyme
    input_fasta = [
        [ id: myid, single_end:false ], // meta map
        file(params.fasta, checkIfExists: true)
    ]

    MAKECMAP_FA2CMAPMULTICOLOR ( input_fasta, params.enzyme )

    MAKECMAP_RENAMECMAPIDS(MAKECMAP_FA2CMAPMULTICOLOR.out.cmap, MAKECMAP_FA2CMAPMULTICOLOR.out.cmapkey)

    MAKECMAP_RENAMECMAPIDS.out.renamedcmap.view()

    //new_thing=MAKECMAP_RENAMECMAPIDS.out.renamedcmap.map{a,b ->b}

    //CMAP2BED(new_thing, params.enzyme, params.sample)

    

}
