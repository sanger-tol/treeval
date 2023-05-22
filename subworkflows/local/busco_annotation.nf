#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of organisms to return
// bigbed files containing alignment data between the input fasta and csv style organism names.
// Input - Assembled genomic fasta file
// Output - A BigBed file per datatype per organism entered via csv style in the yaml.

nextflow.enable.dsl=2

// MODULE IMPORT
include { BUSCO_GENE                        } from './busco_gene'
include { ANCESTRAL                         } from './ancestral_gene'

workflow BUSCO_ANNOTATION {
    take:
    dot_genome          // Channel: [val(meta), [ datafile ]]
    reference_tuple
    assembly_classT
    lineageinfo          // channel: [ meta, lineage_db, lineage ] 
    lineagespath         // channel: [ meta, /path/to/buscoDB, lineage ] 
    buscogene_as         // channel: val(dot_as location)
    ancestral_table      // channel: val(ancestral_table location)
    

    main:
    ch_versions         = Channel.empty()

    //
    // LOGIC: TAKES A SINGLE LIKE CSV STRING AND CONVERTS TO LIST OF VALUES
    //          LIST IS MERGED WITH DATA_DIRECTORY AND ORGANISM_CLASS
    //

    ch_data = reference_tuple
              .combine( lineageinfo )
              .combine( lineagespath)
              .combine( assembly_classT )
              .combine( ancestral_tables)
        
        .branch {
            lep:     it[6]  == 'lep'
            general: it[6]  != 'lep'
        }
        .multiMap {
            refmeta, ref, limeta, li, lpmeta, lp, classT, atable ->
            refT  : tuple(refmeta, ref)
            liT   : tuple(limeta, li)
            lpT   : tuple(lpmeta, lp)
            classT: classT
            atable: atable
            
        }
        .set {ch_busco_data}

    //
    // SUBWORKFLOW: GENERATES GENE ALIGNMENTS FOR RNA, NUCLEAR AND COMPLEMENT_DNA DATA, EMITS BIGBED
    //
    BUSCO_GENE (        ch_busco_data.general.refT,
                        ch_busco_data.general.liT,
                        ch_busco_data.general.lpT,
                        dot_genome,
                        buscogene_as )
    
    ANCESTRAL (         ch_busco_data.lep.refT,
                        ch_busco_data.lep.liT,
                        ch_busco_data.lep.lpT,
                        dot_genome,
                        buscogene_as,
                        ch_busco_data.lep.atable)

    emit:
    buscogene_bigbed        = BUSCO_GENE.out.bigbed
    buscogene_version       = BUSCO_GENE.out.versions
    ancestral_bigbed        = BUSCO_GENE.out.bigbed
    ancestral_version       = BUSCO_GENE.out.versions

}
