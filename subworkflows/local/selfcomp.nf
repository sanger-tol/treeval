#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// MODULE IMPORT
include { MUMMER } from '../../modules/nf-core/modules/mummer/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'
include { SELFCOMP_SPLITFASTA as SPLITFASTA } from '../../modules/sanger-tol/nf-core-modules/selfcomp/splitfasta/main'
include { SELFCOMP_MUMMER2BED } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mummer2bed/main'
include { SELFCOMP_MAPIDS } from '../../modules/sanger-tol/nf-core-modules/selfcomp/mapids/main'
include { CHUNKFASTA } from '../../modules/local/chunkfasta'

workflow SELFCOMP {

    main:
    ch_versions = Channel.empty()

    //def version = '0.001-c3'

    // Inputs
    input_fasta = params.reference
    sample_name = params.assembly.sample
    //take:.genome
    // - assembly.as file
    
    // 1. Generate chromsize value
    //SAMTOOLS_FAIDX([[id:sample_name], input_fasta])
    //ch_fai = SAMTOOLS_FAIDX.out.fai
    //cut -f1,2 ch_fai|sort -k2,2 -nr > my.genome

    //2. Split fasta sequences into 500kb length fasta sequences
    SPLITFASTA([[id:sample_name, single_end: true], input_fasta])
    ch_split_fa = SPLITFASTA.out.fa
    ch_split_agp = SPLITFASTA.out.agp

    ch_split_fa.view()

    // 3. Chunk the fasta file for processing (1Gb chunks)
    CHUNKFASTA(ch_split_fa)
    ch_chunked_fas = CHUNKFASTA.out.fas

    // 4. Run mummer on each chunk, generating .coords files for each.

        // while read p; do
        // j=$(basename $p .fa)
        // echo "/software/grit/bin/mummer3.23 -n -b -c -L -l 400 ${p} split.0.fa split.3.fa split.1.fa split.2.fa split.4.fa split.5.fa > ${j}.mummer" | bsub -q normal -o ${j}.out -e ${j}.err -n 8 -M10000 -R'select[mem>10000] rusage[mem=10000] span[hosts=1]'
        // done < scaffold.fofn

    process CHECK_CHANNELS {

        input:
        tuple val(meta), path(ref), path(query)
        
        output:
        tuple val(meta), path(ref), path(query)
        
        shell:
        """
        echo $meta $ref $query
        """
    }

    //ch_chunked_fas
    //    .combine(ch_chunked_fas)
    //    .view()

    //CHECK_CHANNELS(ch_chunked_fas)

    // Start channel of tuple (meta, path)


    // End channel of tuples

    //MUMMER(ch_chunked_fas)

    // 5. Concatenate mummer files.


    // 6. Run mummer2bed module (.mummer, motiflen) -> bed.
    //MUMMER2BED([[id:sample_name, single_end: true], input_fasta], ch_split_agp)
    //ch_bed= MUMMER2BED.out.bed

    // 7. Run mapids module (.bed, .agp) -> .bed
    //MAPIDS([[id:sample_name, single_end: true], ch_bed])
    //ch_mapped_bed = MAPIDS.out.bed


    // 8. Convert bed to bigbed, requires chromsize and assembly.as
    // /software/grit/tools/ucsc/bedToBigBed -as=new.as -type=bed3+5 -extraIndex=name,qStart,qEnd filtered.bed my.genome rLacAgi1_1_selfcomp.bb
    // Output 
    // selfcomp.bb file


    emit:
    ch_versions

}