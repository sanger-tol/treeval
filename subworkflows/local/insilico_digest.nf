#!/usr/bin/env nextflow
//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

include { MAKECMAP_FA2CMAPMULTICOLOR } from '../../modules/sanger-tol/nf-core-modules/makecmap/fa2cmapmulticolor/main'
include { MAKECMAP_RENAMECMAPIDS } from '../../modules/sanger-tol/nf-core-modules/makecmap/renamecmapids/main'
include { MAKECMAP_CMAP2BED } from '../../modules/sanger-tol/nf-core-modules/makecmap/cmap2bed/main'
include { BED2BIGBED } from '../../modules/local/bed2bigbed'



nextflow.enable.dsl = 2

workflow INSILICO_DIGEST {

    main:
    enzyme = params.enzyme
    sample = params.sample
    sizefile = params.chromsize
    myid = sample + "_" + enzyme
    input_fasta = [
        [ id: myid, single_end:false ], // meta map
        file(params.fasta, checkIfExists: true)
    ]

    MAKECMAP_FA2CMAPMULTICOLOR ( input_fasta, enzyme )

    MAKECMAP_RENAMECMAPIDS(MAKECMAP_FA2CMAPMULTICOLOR.out.cmap, MAKECMAP_FA2CMAPMULTICOLOR.out.cmapkey)

    ch_editedcmap = MAKECMAP_RENAMECMAPIDS.out.renamedcmap

    ch_editedcmap.view()

    MAKECMAP_CMAP2BED( ch_editedcmap, enzyme )

    ch_bed = MAKECMAP_CMAP2BED.out.bedfile

    ch_bed.view()

    BED2BIGBED(ch_bed, sizefile, '-as=../assets/digest.as -type=bed4+1 -extraIndex=length' ${myid}.bigbed)
}
