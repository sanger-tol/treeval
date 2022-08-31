#!/usr/bin/env nextflow
//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

include { MAKECMAP_FA2CMAPMULTICOLOR } from '../../modules/sanger-tol/nf-core-modules/makecmap/fa2cmapmulticolor/main'
include { MAKECMAP_RENAMECMAPIDS } from '../../modules/sanger-tol/nf-core-modules/makecmap/renamecmapids/main'
include { MAKECMAP_CMAP2BED } from '../../modules/sanger-tol/nf-core-modules/makecmap/cmap2bed/main'
include { UCSC_BEDTOBIGBED } from '../../modules/nf-core/modules/ucsc/bedtobigbed/main'



nextflow.enable.dsl = 2

workflow INSILICO_DIGEST {

    main:

    sample = params.sample
    sizefile = params.chromsize
    myid = sample

    ch_enzyme = Channel.of( "bspq1","bsss1","DLE1" )
    ch_versions = Channel.empty()

    input_fasta = [
        [ id: myid, single_end:false ], // meta map
        file(params.fasta, checkIfExists: true)
    ]

    MAKECMAP_FA2CMAPMULTICOLOR ( input_fasta, ch_enzyme )

    ch_cmap    = MAKECMAP_FA2CMAPMULTICOLOR.out.cmap
    ch_cmapkey = MAKECMAP_FA2CMAPMULTICOLOR.out.cmapkey
    ch_version = ch_versions.mix(MAKECMAP_FA2CMAPMULTICOLOR.out.versions)


    ch_cmap_new = ch_cmap
        .map{ meta, cfile  -> tuple([
                                    id  :  cfile.toString().split('_')[-3]
        ], cfile)} 

    ch_cmapkey_new = ch_cmapkey
        .map{ kfile  -> tuple([
                                    id  :  kfile.toString().split('_')[-4]
        ], kfile)}


    ch_join = ch_cmap_new.join(ch_cmapkey_new)
        .map { meta, cfile, kfile -> tuple ([
                                                meta,
                                                cfile
                                                ] ,
                                            kfile)}
 
    MAKECMAP_RENAMECMAPIDS ( ch_join.map { it[0] }, ch_join.map { it[1] } )
    ch_version = ch_versions.mix(MAKECMAP_RENAMECMAPIDS.out.versions)

    ch_renamedcmap = MAKECMAP_RENAMECMAPIDS.out.renamedcmap

    MAKECMAP_CMAP2BED ( ch_renamedcmap, ch_renamedcmap.map { it[0].id } )
    ch_version = ch_versions.mix(MAKECMAP_CMAP2BED.out.versions)

    ch_bedfile = MAKECMAP_CMAP2BED.out.bedfile

    UCSC_BEDTOBIGBED ( ch_bedfile, sizefile)
    ch_version = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    versions = ch_version

}
