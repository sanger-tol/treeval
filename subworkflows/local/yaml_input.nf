#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

workflow YAML_INPUT {
    take:
    input_file      // params.input
    workflow_name   // params.entry

    main:
    ch_versions = Channel.empty()

    input_file
        .map { file -> readYAML(file) }
        .set { yamlfile }

    Channel.of( workflow_name )
        .set{ workflow_id }

    //
    // LOGIC: PARSES THE TOP LEVEL OF YAML VALUES
    //
    yamlfile
        .flatten()
        .combine( workflow_id )
        .multiMap { data, id ->
                assembly:               ( data.assembly )
                assembly_reads:         ( data.assem_reads )
                hic_data:               ( data.hic_data )
                kmer_profile:           ( data.kmer_profile )
                reference:              ( file(data.reference_file, checkIfExists: true) )
                alignment:              ( id == "FULL" ? data.alignment : "" )
                self_comp:              ( id == "FULL" ? data.self_comp : "" )
                synteny:                ( id == "FULL" ? data.synteny   : "" )
                intron:                 ( id == "FULL" ? data.intron    : "" )
                busco_gene:             ( data.busco )
                teloseq:                ( data.telomere )
                map_order:              ( data.map_order)
        }
        .set{ group }

    group.synteny.view{"INPUT: $it"}

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
        .self_comp
        .combine( workflow_id )
        .multiMap { data, id ->
                    motif_len:          (id == "FULL" ? data.motif_len          : "")
                    mummer_chunk:       (id == "FULL" ? data.mummer_chunk       : "")
        }
        .set{ selfcomp_data }

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

    if ( assem_reads.read_type.filter { it == "hifi" } || assem_reads.read_type.filter { it == "clr" } || assem_reads.read_type.filter { it == "ont" } ) {
        tolid_version
            .combine( assem_reads.read_type )
            .combine( assem_reads.read_data )
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
            .combine( assem_reads.read_data )
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
        .combine( hic.hic_cram )
        .combine( hic.hic_aligner )
        .map { sample, data, aligner ->
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
    assembly_id                      = tolid_version
    reference_ch                     = ref_ch
    map_order_ch                     = group.map_order

    read_ch                          = read_ch

    kmer_prof_file                   = kmer_prof

    hic_reads_ch                     = hic_ch
    supp_reads_ch                    = supplement_ch

    // align_data_dir                   = alignment_data.data_dir
    // align_geneset                    = alignment_data.geneset_id
    // align_common                     = alignment_data.common_name
    align_genesets                    = alignment_data.genesets

    motif_len                        = selfcomp_data.motif_len
    mummer_chunk                     = selfcomp_data.mummer_chunk

    synteny_paths                    = group.synteny

    intron_size                      = intron_size.size

    teloseq                          = teloseq.teloseq

    lineageinfo                      = busco_lineage.lineage
    lineagespath                     = busco_lineage.lineages_path

    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
