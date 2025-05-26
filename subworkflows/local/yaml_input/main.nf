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
                    assem_level:        data.assem_level        // String
                    assem_version:      data.assem_version      // Number
                    sample_id:          data.sample_id          // String
                    latin_name:         data.latin_name         // String
                    defined_class:      data.defined_class      // String
                    project_id:         data.project_id         // String
            }
        .set { assembly_data }


    group
        .assembly_reads
        .multiMap { data ->
                    read_type:          data.read_type          // String
                    read_data:          file(data.read_data, checkIfExists: true, type: 'dir')
                    supplement:         data.supplementary_data // ?
        }
        .set { assem_reads }


    group
        .hic_data
        .multiMap { data ->
                    hic_cram:          file(data.hic_cram, checkIfExists: true, type: 'dir')
                    hic_aligner:       data.hic_aligner          // String
        }
        .set { hic }


    group
        .kmer_profile
        .multiMap { data ->
                    length:             data.kmer_length         // Number
                    dir:                file(data.dir, checkIfExists: true, type: 'dir')
        }
        .set { kmer_profiling }


    group
        .alignment
        .combine(workflow_id)
        .multiMap{  data, id ->
                    genesets:           (id == "FULL" || id == "JBROWSE" ? data.genesets.collect{ geneset_path -> file(geneset_path, checkIfExists: true) } : "")
        }
        .set{ alignment_data }


    group
        .intron
        .combine( workflow_id )
        .multiMap { data, id ->
                    size:			    (id == "FULL" ? data.size               : "") // String
        }
        .set { intron_size }


    group
        .teloseq
        .multiMap { data ->
                    teloseq:            data.teloseq               // String
        }
        .set { teloseq }


    group
        .busco_gene
        .multiMap { data ->
                    lineage:			data.lineage                // String
                    lineages_path:		file(data.lineages_path, checkIfExists: true, type:'dir')
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

    tolid_version
        .combine( assem_reads.read_type )
        .combine( assem_reads.read_data )
        .filter{ _sample, type, _data -> type in ['hifi', 'clr', 'ont', 'illumina'] }
        .map{ sample, type, data ->
            tuple(  [   id              : sample,
                        single_end      : type != "illumina",
                        read_type       : type
                    ],
                    data
            )
        }
    .set { read_ch }


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

    align_genesets                   = alignment_data.genesets

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
