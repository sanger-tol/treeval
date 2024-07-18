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
    sizefile        // Channel: tuple [ val(meta), path(my.genome_file) ]
    reference       // Channel: tuple [ val(meta), path(reference_file) ]
    ch_enzyme       // Channel: val( "bspq1","bsss1","DLE1" )
    dot_as          // Channel: val(dot_as location)

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: COMBINES REFERENCE TUPLE WITH ENZYME CHANNEL
    //        MULTIMAP INTO TWO CHANNELS SO THERE IS REFERENCE * ENZYME CHANNELS
    //
    reference
        .map {meta, data ->
            tuple(
                [   id          : meta.id,
                    single_end  : false     ],
                file( data )
            )
    }
    .set {input_fasta}

    input_fasta
        .combine(ch_enzyme)
        .multiMap {meta, reference, enzyme_id ->
            fasta       : tuple(    meta,
                                    reference
                            )
            enzyme      : enzyme_id
            }
        .set {fa2c_input}

    //
    // MODULE: CONVERTS FASTA INTO A COLOUR-AWARE BIONANO CMAP FORMAT
    //         EMITS FILES CONTAINING INDEX_IDs AND ORIGINAL_GENOMIC_LOCATIONS
    //
    MAKECMAP_FA2CMAPMULTICOLOR (
        fa2c_input.fasta,
        fa2c_input.enzyme
    )
    ch_versions         = ch_versions.mix(MAKECMAP_FA2CMAPMULTICOLOR.out.versions)

    //
    // LOGIC: CREATES A TUPLE CONTAINING THE CMAP AND ORIGINAL GENOMIC LOCATIONS
    //
    MAKECMAP_FA2CMAPMULTICOLOR.out.cmap
        .map{ meta, cfile  ->
            tuple(
                [id    :  cfile.toString().split('_')[-3]],
                cfile
            )
        }
        .set { ch_cmap_new }

    MAKECMAP_FA2CMAPMULTICOLOR.out.cmapkey
        .map{ kfile  ->
            tuple(
                [id    :  kfile.toString().split('_')[-4]],
                kfile
            )
        }
        .set {ch_cmapkey_new}


    ch_cmap_new
        .join(ch_cmapkey_new)
        .multiMap {meta, cfile, kfile ->
            cmap        : tuple( meta, cfile)
            key_file    : kfile
        }

        .set {ch_join}

    //
    // MODULE: RENAME CMAP IDs FROM BIONANO IDX TO ORIGINAL GENOMIC LOCATIONS
    //         EMITS RENAMED CMAP
    //
    MAKECMAP_RENAMECMAPIDS (
        ch_join.cmap,
        ch_join.key_file
    )
    ch_versions         = ch_versions.mix(MAKECMAP_RENAMECMAPIDS.out.versions)

    MAKECMAP_RENAMECMAPIDS.out.renamedcmap
        .multiMap {meta, file ->
            full        : tuple ( meta, file )
            sample      : meta.id
        }
        .set {ch_renamedcmap}

    //
    // MODULE: CONVERT CMAP FILE INTO BED FILE
    //         EMITS BED FILE
    //
    MAKECMAP_CMAP2BED (
        ch_renamedcmap.full,
        ch_renamedcmap.sample
    )
    ch_versions         = ch_versions.mix(MAKECMAP_CMAP2BED.out.versions)

    MAKECMAP_CMAP2BED.out.bedfile
        .combine(sizefile)
        .combine(dot_as)
        .multiMap {meta, bed, meta_2, dot_genome, as_file ->
            bed_tuple   : tuple(meta, bed)
            genome_file : dot_genome
            autosql     : as_file
        }
        .set {combined_ch}

    //
    // MODULE: CONVERT ABOVE BED INTO BIGBED WITH ADDITIONAL AS FILE
    //         EMITS BIGBED FILE
    //
    UCSC_BEDTOBIGBED (
        combined_ch.bed_tuple,
        combined_ch.genome_file,
        combined_ch.autosql
    )
    ch_versions         = ch_versions.mix(UCSC_BEDTOBIGBED.out.versions)

    emit:
    insilico_digest_bb  = UCSC_BEDTOBIGBED.out.bigbed
    versions            = ch_versions.ifEmpty(null)
}
