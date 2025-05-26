#!/usr/bin/env nextflow

workflow YAML_INPUT {
    take:
    input_file      // params.input
    workflow_name   // params.mode

    main:
    ch_versions = Channel.empty()

    yamlfile    =   Channel.of(
                        readYAML(input_file)
                    )

    workflow_id =   Channel.of(
                        workflow_name
                    )


    //
    // LOGIC: PARSES THE TOP LEVEL OF YAML VALUES
    //

    yamlfile
        .flatten()
        .multiMap { data ->
                assembly:               ( data.assembly )
                assembly_reads:         ( data.assem_reads )
                hic_data:               ( data.hic_data )
                kmer_profile:           ( data.kmer_profile )
                reference:              ( file(data.reference_file, checkIfExists: true) )
                alignment:              ( params.mode == "FULL" ? data.alignment : "" )
                synteny:                ( data.synteny ? data.synteny   : "" )
                intron:                 ( params.mode == "FULL" ? data.intron    : "" )
                busco_gene:             ( data.busco )
                teloseq:                ( data.telomere )
                map_order:              ( data.map_order)
        }
        .set{ group }


    //
    // LOGIC: PARSES THE SECOND LEVEL OF YAML VALUES PER ABOVE OUTPUT CHANNEL
    //
    group
        .assembly
        .multiMap { data ->
                    assem_level:        data.assem_level
                    assem_version:      data.assem_version
                    sample_id:          data.sample_id
                    latin_name:         data.latin_name
                    defined_class:      data.defined_class
                    project_id:         data.project_id
            }
        .set { assembly_data }


    group
        .assembly_reads
        .multiMap { data ->
                    read_type:          data.read_type
                    read_data:          data.read_data
                    supplement:         data.supplementary_data
        }
        .set { assem_reads }


    group
        .hic_data
        .multiMap { data ->
                    hic_cram:          data.hic_cram
                    hic_aligner:       data.hic_aligner
        }
        .set { hic }


    group
        .kmer_profile
        .multiMap { data ->
                    length:             data.kmer_length
                    dir:                data.dir
        }
        .set { kmer_profiling }


    group
        .alignment
        .combine(workflow_id)
        .multiMap{  data, id ->
                    genesets:           (id == "FULL" || id == "JBROWSE" ? data.genesets           : "")
        }
        .set{ alignment_data }


    group
        .intron
        .combine( workflow_id )
        .multiMap { data, id ->
                    size:			    (id == "FULL" ? data.size               : "")
        }
        .set { intron_size }


    group
        .teloseq
        .multiMap { data ->
                    teloseq:            data.teloseq
        }
        .set { teloseq }


    group
        .busco_gene
        .multiMap { data ->
                    lineage:			data.lineage
                    lineages_path:		data.lineages_path
        }
        .set { busco_lineage }


    //
    // LOGIC: COMBINE SOME CHANNELS INTO VALUES REQUIRED DOWNSTREAM
    //
    assembly_data.sample_id
        .combine( assembly_data.assem_version )
        .map { it1, it2 ->
            ("${it1}_${it2}")}
        .set { tolid_version }


    tolid_version
        .combine( group.reference )
        .combine( assembly_data.defined_class )
        .combine( assembly_data.project_id )
        .map { sample, ref_file, defined_class, project ->
            tuple(  [   id:             sample,
                        class:          defined_class,
                        project_type:   project
                    ],
                    ref_file
            )
        }
        .set { ref_ch }


    //
    // LOGIC: COLLECT ARRAY INTO LIST
    //
    ch_collected_reads  = assem_reads.read_data
                            .collect()
                            .map{ files -> [files] }

    ch_collected_hic    = hic.hic_cram
                            .collect()
                            .map{ files -> [files] }

    //
    // LOGIC: FIND THE FOFN IF THERE IS ONE AND PULL OUT THE FILE NAME PER LINE
    //
    filtered_pb_ch      = ch_collected_reads
        .flatten()
        .transpose()
        .filter{it.toString().contains(".fofn")}
        .map{it -> file(it).text.split('\n').collect { it.trim() }}

    filtered_hic_ch     = ch_collected_hic
        .flatten()
        .transpose()
        .filter{it.toString().contains(".fofn")}
        .map{it -> file(it).text.split('\n').collect { it.trim() }}

    //
    // LOGIC: FIND ERRONEOUS LINES IN THE INPUT FILES (NOT YET THE PROCESSED FILES) COUNT ERRORS AND ERROR IF GREATER THAN 0
    //
    ch_collected_reads
        .flatten()
        .transpose()
        .filter{ !it.toString().contains(".fofn") && !it.toString().contains(".fasta.gz") && !it.toString().contains("fa.gz") }
        .count().map { n ->
            if (n > 0) {
                exit 1, "One of the input longread files does not match fa.gz, fasta.gz, fofn."
            }
        }

    ch_collected_hic
        .flatten()
        .transpose()
        .filter{ !it.toString().contains(".fofn") && !it.toString().contains(".cram")}
        .count().map { n ->
            if (n > 0) {
                exit 1, "One of the input hic files does not match cram or fofn."
            }
        }


    //
    // LOGIC: IF THERE IS A FOFN, MERGE WITH INPUT (WHICH MAY CONTAIN OTHER READS) AND THEN COLLECT AS LIST.
    //
    if (filtered_pb_ch) {
        ch_collected_reads
            .flatten()
            .transpose()
            .filter{!it.toString().contains(".fofn")}
            .concat(filtered_pb_ch.flatten())
            .collect()
            .map{ files -> [files]}
            .set {final_pb_reads}

        // NOTE: This will check the processed files used as input (pooled from fofn contents AND input array)
        final_pb_reads.flatten().transpose().filter{!it.toString().contains(".fasta.gz")}.count().map { n ->
            if (n > 0) {
                exit 1, "One of the input read files does not match `fa.gz`, `fasta.gz`. CHECK YOUR FOFN CONTENTS"
            }
        }
    } else {
        // NOTE: IF NO FOFN JUST OUTPUT THE INPUT CHANNEL
        final_pb_reads = ch_collected_reads
    }

    if (filtered_hic_ch) {
        ch_collected_hic
            .flatten()
            .transpose()
            .filter{!it.toString().contains(".fofn")}
            .concat(filtered_hic_ch.flatten())
            .collect()
            .map{ files -> [files]}
            .set {final_hic_reads}

        // NOTE: This will check the processed files used as input (pooled from fofn contents AND input array)
        final_hic_reads.flatten().transpose().filter{!it.toString().contains(".cram")}.count().map { n ->
            if (n > 0) {
                exit 1, "One of the input read files does not match `cram`. CHECK YOUR FOFN CONTENTS"
            }
        }
    } else {
        // NOTE: IF NO FOFN JUST OUTPUT THE INPUT CHANNEL
        final_hic_reads = ch_collected_hic
    }


    if ( assem_reads.read_type.filter { it == "hifi" } || assem_reads.read_type.filter { it == "clr" } || assem_reads.read_type.filter { it == "ont" } ) {
        tolid_version
            .combine( assem_reads.read_type )
            .combine( final_pb_reads )
            .map{ sample, type, data ->
                tuple(  [   id              : sample,
                            single_end      : true,
                            read_type       : type
                        ],
                        data
                )
            }
        .set { read_ch }
    }
    else if ( assem_reads.read_type.filter { it == "illumina" } ) {
        tolid_version
            .combine( assem_reads.read_type )
            .combine( final_pb_reads )
            .map{ sample, type, data ->
                tuple(  [   id              : sample,
                            single_end      : false,
                            read_type       : type
                        ],
                        data
                )
            }
        .set { read_ch }
    }


    tolid_version
        .combine( hic.hic_aligner )
        .combine( final_hic_reads )
        .map { sample, aligner, data ->
            tuple(  [   id: sample,
                        aligner: aligner  ],
                    data
            )
        }
        .set { hic_ch }


    tolid_version
        .combine( assem_reads.supplement )
        .map { sample, data ->
            tuple(  [   id: sample  ],
                    data
            )
        }
        .set { supplement_ch }


    tolid_version
        .combine ( assembly_data.sample_id )
        .combine ( kmer_profiling.length )
        .combine ( kmer_profiling.dir )
        .map { sample, sample_id, kmer_len, dir ->
            tuple(  [   id: sample,
                        kmer: kmer_len  ],
                    file("${dir}/k${kmer_len}/${sample_id}.k${kmer_len}.ktab") // Don't check for existence yet
            )
        }
        .set { kmer_prof }

    emit:
    ch_assembly_id      = tolid_version
    ch_reference        = ref_ch
    ch_map_order        = group.map_order

    ch_assem_reads      = read_ch

    ch_kmer_prof_file   = kmer_prof

    ch_hic_reads        = hic_ch
    ch_supp_reads       = supplement_ch

    ch_align_genesets   = alignment_data.genesets

    ch_synteny_paths    = group.synteny

    ch_intron_size      = intron_size.size

    ch_teloseq          = teloseq.teloseq

    ch_lineageinfo      = busco_lineage.lineage
    ch_lineagespath     = busco_lineage.lineages_path

    versions            = ch_versions
}

def readYAML( yamlfile ) {
    return new org.yaml.snakeyaml.Yaml().load( new FileReader( yamlfile.toString() ) )
}
