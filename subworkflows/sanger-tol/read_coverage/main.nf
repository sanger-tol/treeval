include { MINIMAP2_ALIGN     } from '../../../modules/sanger-tol/minimap2/align/main'
include { FIND_CONCATENATE   } from '../../../modules/sanger-tol/find/concatenate/main'
include { BEDTOOLS_GENOMECOV } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow READ_COVERAGE {

    take:
    ch_reads      // channel: [ val(meta), [ path(reads) ] ]
    ch_reference  // channel: [ val(meta2), path(fasta) ]
    ch_chromsizes // channel: [ path(chromsizes) ]

    main:

    //
    // MODULE: Normalise reads and pair each read with its matching reference (by meta id)
    //
    ch_reads
        .map { meta, reads ->
            def read_list = (reads instanceof List) ? reads : [reads]
            tuple(meta, read_list)
        }
        .flatMap { meta, read_list ->
            read_list.withIndex().collect { read_file, idx ->
                def meta_with_read_idx = meta + [read_idx: idx]
                tuple(meta_with_read_idx, read_file)
            }
        }
        .combine(ch_reference)
        .multiMap { meta, read_file, meta2, reference ->
            reads: tuple(meta, read_file)
            reference: tuple(meta2, reference)
        }
        .set { ch_align_split }


    //
    // MODULE: Run minimap2 and emit one BED per read file
    //
    MINIMAP2_ALIGN (
        ch_align_split.reads,
        ch_align_split.reference,
        false,
        [],
        false,
        false,
        true
    )

    // NOTE: Group per-sample PAF/BED outputs into lists for downstream concatenation
    ch_paf_bed_grouped = MINIMAP2_ALIGN.out.bed
        .groupTuple(by: 0)
        .map { meta, paf_bed_files -> [ meta, paf_bed_files.sort { f -> f.name } ] }


    //
    // MODULE: Merge all per-read BED outputs into one BED per sample and sort it
    //
    FIND_CONCATENATE(ch_paf_bed_grouped, true)

    find_concat_out = FIND_CONCATENATE.out.file_out
        .map{ meta, bed -> [[id: meta.id], bed, 1] }


    //
    // MODULE: Generate BedGraph from merged, sorted BED
    //
    BEDTOOLS_GENOMECOV (
        find_concat_out,
        ch_chromsizes,
        'bedgraph',
        true
    )


    //
    // MODULE: Convert to BigWig (Compulsory)
    //
    UCSC_BEDGRAPHTOBIGWIG (
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_chromsizes
    )


    emit:
    bigwig         = UCSC_BEDGRAPHTOBIGWIG.out.bigwig
    bedgraph       = BEDTOOLS_GENOMECOV.out.genomecov
}
