#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

workflow YAML_INPUT {
    take:
    input_file  // params.input

    main:
    ch_versions = Channel.empty()

    input_file
        .map { file -> readYAML(file) }
        .set { yamlfile }

    //
    // LOGIC: PARSES THE TOP LEVEL OF YAML VALUES
    //
    yamlfile
        .flatten()
        .multiMap { data ->
                assembly:               ( data.assembly )
                assembly_reads:         ( data.assem_reads )
                reference:              ( file(data.reference_file, checkIfExists: true) )
                alignment:              ( data.alignment )
                self_comp:              ( data.self_comp )
                synteny:                ( data.synteny )
                intron:                 ( data.intron )
                busco_gene:             ( data.busco )
                teloseq:                ( data.telomere )
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
                    longread_type:      data.longread_type
                    longread_data:      data.longread_data
                    hic:                data.hic_data
                    supplement:         data.supplementary_data
        }
        .set { assem_reads }

    group
        .alignment
        .multiMap { data ->
                    data_dir:           data.data_dir
                    common_name:        data.common_name
                    geneset_id:         data.geneset_id
        }
        .set{ alignment_data }

    group
        .self_comp
        .multiMap { data ->
                    motif_len:          data.motif_len
                    mummer_chunk:       data.mummer_chunk
        }
        .set{ selfcomp_data }

    group
        .synteny
        .multiMap { data ->
                    synteny_genome:     data.synteny_genome_path
        }
        .set{ synteny_data }

    group
        .intron
        .multiMap { data ->
                    size:			    data.size
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
        .set { tolid_version}

    assembly_data.sample_id
        .combine( assembly_data.assem_version )
        .combine( group.reference )
        .combine( assembly_data.defined_class )
        .combine( assembly_data.project_id )
        .map { sample, version, ref_file, defined_class, project ->
            tuple(  [   id:             "${sample}_${version}",
                        class:          defined_class,
                        project_type:   project
                    ],
                    ref_file
            )
        }
        .set { ref_ch }

    assembly_data.sample_id
        .combine( assem_reads.longread_type )
        .combine( assem_reads.longread_data )
        .map{ sample, type, data ->
            tuple(  [   id              : sample,
                        single_end      : true,
                        longread_type   : type
                    ],
                    data
            )
        }
        .set { longread_ch }

    assembly_data.sample_id
        .combine( assem_reads.hic )
        .map { sample, data ->
            tuple(  [   id: sample  ],
                    data
            )
        }
        .set { hic_ch }

    assembly_data.sample_id
        .combine( assem_reads.supplement )
        .map { sample, data ->
            tuple(  [   id: sample  ],
                    data
            )
        }
        .set { supplement_ch }

    emit:
    reference                        = ref_ch
    assembly_level                   = assembly_data.assem_level

    assembly_id                      = tolid_version
    assembly_classT                  = assembly_data.defined_class
    assembly_level                   = assembly_data.assem_level
    assembly_asmVer                  = assembly_data.assem_version
    assembly_ttype                   = assembly_data.project_id

    longreads_new                        = longread_ch
    hic_reads_new                        = hic_ch
    supp_reads_new                       = supplement_ch

    pacbio_reads                     = assem_reads.longread_data
    hic_reads                        = assem_reads.hic
    supp_reads                       = assem_reads.supplement

    align_data_dir                   = alignment_data.data_dir
    align_geneset                    = alignment_data.geneset_id
    align_common                     = alignment_data.common_name

    motif_len                        = selfcomp_data.motif_len
    mummer_chunk                     = selfcomp_data.mummer_chunk

    synteny_path                     = synteny_data.synteny_genome

    intron_size                      = intron_size.size

    teloseq                          = teloseq.teloseq

    lineageinfo                      = busco_lineage.lineage
    lineagespath                     = busco_lineage.lineages_path

    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
