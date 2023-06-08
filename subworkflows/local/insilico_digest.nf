#!/usr/bin/env nextflow
//
// The subworkflow takes an assembly fasta file and produce binano insilico digest cut sites track in bigbed
// Input - genome fasta
// Output - bigbed

//
// MODULE IMPORT BLOCK
//
include { MAKECMAP_FA2CMAPMULTICOLOR    } from '../../modules/local/makecmap_fa2cmapmulticolor'
include { MAKECMAP_RENAMECMAPIDS        } from '../../modules/local/makecmap_renamecmapids'
include { MAKECMAP_CMAP2BED             } from '../../modules/local/makecmap_cmap2bed'
include { UCSC_BEDTOBIGBED              } from '../../modules/nf-core/ucsc/bedtobigbed/main'

workflow INSILICO_DIGEST {
    take:
    myid            // Channel val(sample_id)
    sizefile        // Channel [ val(meta), path(my.genome_file) ]
    sample          // Channel [ val(meta), path(reference_file) ]
    ch_enzyme       // Channel val( "bspq1","bsss1","DLE1" )
    dot_as          // Channel val(dot_as location)

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: COMBINES REFERENCE TUPLE WITH ENZYME CHANNEL
    //        MULTIMAP INTO TWO CHANNELS SO THERE IS REFERENCE * ENZYME CHANNELS
    //
    input_fasta = sample.map { data -> 
                                tuple([
                                    id               : data[0].id,
                                    single_end       : false
                                    ],
                                    file(data[1])
                                )}

    input_fasta
        .combine(ch_enzyme)
        .multiMap { data -> 
            fasta:      tuple( data[0],
                                data[1]
                            )
            enzyme:     data[2]
            }
        .set { fa2c_input } 

    //
    // MODULE: CONVERTS FASTA INTO A COLOUR-AWARE BIONANO CMAP FORMAT
    //         EMITS FILES CONTAINING INDEX_IDs AND ORIGINAL_GENOMIC_LOCATIONS
    //
    MAKECMAP_FA2CMAPMULTICOLOR ( fa2c_input.fasta, fa2c_input.enzyme )

    ch_cmap             = MAKECMAP_FA2CMAPMULTICOLOR.out.cmap
    ch_cmapkey          = MAKECMAP_FA2CMAPMULTICOLOR.out.cmapkey
    ch_versions         = ch_versions.mix(MAKECMAP_FA2CMAPMULTICOLOR.out.versions)

    //
    // LOGIC: CREATES A TUPLE CONTAINING THE CMAP AND ORIGINAL GENOMIC LOCATIONS
    //
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
 
    //
    // MODULE: RENAME CMAP IDs FROM BIONANO IDX TO ORIGINAL GENOMIC LOCATIONS
    //         EMITS RENAMED CMAP
    //
    MAKECMAP_RENAMECMAPIDS ( ch_join.map { it[0] }, ch_join.map { it[1] } )
    ch_versions         = ch_versions.mix(MAKECMAP_RENAMECMAPIDS.out.versions)

    ch_renamedcmap      = MAKECMAP_RENAMECMAPIDS.out.renamedcmap

    //
    // MODULE: CONVERT CMAP FILE INTO BED FILE
    //         EMITS BED FILE
    //
    MAKECMAP_CMAP2BED ( ch_renamedcmap, ch_renamedcmap.map { it[0].id } )
    ch_versions         = ch_versions.mix(MAKECMAP_CMAP2BED.out.versions)

    ch_bedfile          = MAKECMAP_CMAP2BED.out.bedfile
    combined_ch         = ch_bedfile
                            .combine(sizefile)
                            .combine(dot_as)
    
    //
    // MODULE: CONVERT ABOVE BED INTO BIGBED WITH ADDITIONAL AS FILE
    //         EMITS BIGBED FILE
    //
    UCSC_BEDTOBIGBED (  combined_ch.map { [it[0], it[1]] },
                        combined_ch.map { it[3] },
                        combined_ch.map { it[4] })
    ch_versions         = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    insilico_digest_bb  = UCSC_BEDTOBIGBED.out.bigbed
    versions            = ch_versions.ifEmpty(null)
}
