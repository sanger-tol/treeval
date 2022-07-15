//
// Check for synteny by aligning to fasta to reference genomes.
//

include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'

workflow SYNTENY {
    take:
    input_fasta // channel: /path/to/fasta
    clade // channel [val(meta), [ val ]]

    main:
    ch_versions = Channel.empty()

    //TO MOVE TO CONFIG
    ch_fasta = Channel.fromPath('/nfs/team135/dp24/treeval_testdata/synteny_data/birds/bTaeGut2.fasta')
    genome_path = '/nfs/team135/dp24/treeval_testdata/synteny_data'

    // If no reference genomes returned for clade then no .PAF outputted.
    clade_val = clade 
    ch_genomes = Channel.fromPath("$genome_path/birds/*.fasta")
    
    // MINIMAP2_ALIGN 
    //  Input { [ meta, reads ], reference, bam_format, cigar_paf_format, cigar_bam }
    //  Output { meta, paf, bam, versions}

    MINIMAP2_ALIGN ( [[id:'bTaeGut2', single_end:false ], ch_fasta], ch_genomes, false, true, false )
    ch_align_paf = MINIMAP2_ALIGN.out.paf.first()

    //publishdir
    ch_results = SYNTENY.out.ch_align_paf
    //ch_results = SYNTENY.out.ch_align_paf
    ch_results.collectFile(name:'paftest.txt')
    //ch_versions = ch_versions.mix(SYNTENY.out.versions)
    //ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    emit:
    ch_align_paf
    //genomes = ch_genomes
    //clade_val
    //ch_versions


}