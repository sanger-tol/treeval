import org.yaml.snakeyaml.Yaml

workflow INPUT_READ {
    take:
    input_file

    main:
    ch_versions = Channel.empty()

    input_file
        .map { file -> readYAML(file) }
        .set { yamlfile }

    // Parse top layer of yaml
    yamlfile
        .flatten()
        .multiMap { data -> 
                assembly:               ( data.assembly )
                assembly_reads:         ( data.assem_reads )
                reference:              ( file(data.reference_file) )
                alignment:              ( data.alignment )
                self_comp:              ( data.self_comp )
                synteny:                ( data.synteny )
                intron:                 data.intron_size
        }
        .set{ group }

    // Parse 2nd layer
    group
        .assembly
        .multiMap { data ->
                    level:              data.level
                    size_c:             data.sizeClass
                    sample_id:          data.sample_id
                    classT:             data.classT
                    asmVersion:         data.asmVersion
                    dbVersion:          data.dbVersion
                    gevalType:          data.gevalType
            }
        .set { assembly_data }

    group
        .assem_reads
        .multiMap { data -> 
            pacbio:                     data.pacbio
            hic:                        data.hic
            supplement:                 data.supplementary
        }
        .set { assembly_reads }

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
                motif_len:              data.motif_len
                mummer_chunk:           data.mummer_chunk
        }
        .set{ selfcomp_data }

    group
        .synteny
        .multiMap { data -> 
                synteny_genome:         data.synteny_genome_path
        }
        .set{ synteny_data }

    emit:
    assembly_id                      = assembly_data.sample_id
    assembly_sizeClass               = assembly_data.size_c
    assembly_classT                  = assembly_data.classT
    assembly_level                   = assembly_data.level
    assembly_asmVer                  = assembly_data.asmVersion
    assembly_dbVer                   = assembly_data.dbVersion
    assembly_gtype                   = assembly_data.gevalType

    pacbio_reads                     = assembly_reads.pacbio
    hic_reads                        = assembly_reads.hic
    supp_reads                       = assembly_reads.supplement

    reference                        = group.reference

    align_data_dir                   = alignment_data.data_dir
    align_geneset                    = alignment_data.geneset
    align_common                     = alignment_data.geneset

    motif_len                        = selfcomp_data.motif_len
    mummer_chunk                     = selfcomp_data.mummer_chunk

    synteny_path                     = synteny_data.synteny_genome

    intron_size                      = group.intron

    versions                         = ch_versions.ifEmpty(null)
}

def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}