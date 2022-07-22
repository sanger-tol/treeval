//
// Check for synteny by aligning to fasta to reference genomes.
//
include { MINIMAP2_ALIGN } from '../../modules/nf-core/modules/minimap2/align/main'

workflow SYNTENY {

    main:
    ch_versions = Channel.empty()

    // If no reference genomes returned for class then no .PAF outputted.
    ch_genomes = Channel.fromPath("$params.synteny_genome_path/$params.assembly.class/*.fasta")

    MINIMAP2_ALIGN([[id:params.assembly.sample, single_end: true], params.reference], ch_genomes, false, true, false)

    ch_paf = MINIMAP2_ALIGN.out.paf
    ch_versions = MINIMAP2_ALIGN.out.versions
    
    process OUTPUT_SYNTENY_FILES {
        publishDir params.outdir

        input:
        tuple val(meta), path(paf_file)
        
        output:
        path paf_file
        
        shell:
        """
        cp $paf_file paf_file
        """
    }
    OUTPUT_SYNTENY_FILES(ch_paf)
    
    emit:
    ch_paf
    ch_versions
}