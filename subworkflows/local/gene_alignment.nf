//
// Check input samplesheet and get read channels
//

// MODULES
include { BLAST_MAKEBLASTDB } from './modules/nf-core/modules/blast/makeblastdb/main'

// INPUT
// params.fasta
// params.tolid
// params.alignment_data

workflow GENE_ALIGNMENT {
    take:
    input_fasta         // channel: [ val(meta), [fasta]]
    alignment_folders   // channel: path(val)

    main:
    ch_versions = channel.empty()

    //
    // Something to check alignment data folders
    // Python script to return? 

    //
    // Makeblastdb - creating a blast db of the input fasta
    //
    BLAST_MAKEBLASTDB ( input_fasta )
    
    //
    // blastn or blastx depending on datatype
    //
    

    //
    // collect output from blast, cat and filter on 90% match
    //


    //
    // Generate .genome file with size of all chromosomes
    //

    //
    // Pull fresh assembly.as file for datatype of alignment data 
    //

    //
    // make bigbed
    //

}