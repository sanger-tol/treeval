#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

workflow YAML_INPUT {
    take:
    input_file  // input_yaml_from_commandline

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
                reference:              ( file(data.reference_file) )
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
                    level:              data.level
                    sample_id:          data.sample_id
                    classT:             data.classT
                    asmVersion:         data.asmVersion
                    dbVersion:          data.dbVersion
                    gevalType:          data.gevalType
            }
        .set { assembly_data }

    group
        .assembly_reads
        .multiMap { data ->
                    pacbio:             data.pacbio
                    hic:                data.hic
                    supplement:         data.supplementary
        }
        .set { assem_reads }

    group
        .alignment
        .multiMap { data ->
                    data_dir:           data.data_dir
                    common_name:        data.common_name
                    geneset:            data.geneset
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

    assembly_data.sample_id
        .combine( assembly_data.asmVersion )
        .map { it1, it2 ->
            ("${it1}_${it2}")}
        .set { tolid_version}

    emit:
    assembly_id                      = tolid_version
    assembly_classT                  = assembly_data.classT
    assembly_level                   = assembly_data.level
    assembly_asmVer                  = assembly_data.asmVersion
    assembly_dbVer                   = assembly_data.dbVersion
    assembly_ttype                   = assembly_data.gevalType

    pacbio_reads                     = assem_reads.pacbio
    hic_reads                        = assem_reads.hic
    supp_reads                       = assem_reads.supplement

    reference                        = group.reference

    align_data_dir                   = alignment_data.data_dir
    align_geneset                    = alignment_data.geneset
    align_common                     = alignment_data.geneset

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
