#!/usr/bin/env nextflow

workflow YAML_INPUT {
    take:
    input_file    // String: params.input
    workflow_name // String: params.mode

    main:
    ch_versions = Channel.empty()

    Channel.value(input_file)
        .map { file -> readYAML(file) }
        .flatten()
        .multiMap { data ->
            def id = workflow_name
            def tolid_ver = "${data.assembly.sample_id}_${data.assembly.assem_version}"
            def kmer_len = data.kmer_profile.kmer_length

            // emit:
            tolid_version: tolid_ver
            reference: tuple(
                [
                    id: tolid_ver,
                    class: data.assembly.defined_class,
                    project_type: data.assembly.project_id,
                ],
                file(data.reference_file, checkIfExists: true),
            )
            map_order: data.map_order
            read_ch: GET_VALIDATED_CHANNEL(
                            "longread",
                            tolid_ver,
                            data.assem_reads.read_type,
                            data.assembly.defined_class,
                            data.assembly.project_id,
                            data.assem_reads.read_data.collect()
                        )
            kmer_prof: tuple(
                [
                    id: tolid_ver,
                    kmer: kmer_len,
                ],
                file("${data.kmer_profile.dir}/k${kmer_len}/${data.assembly.sample_id}.k${kmer_len}.ktab"),
            )
            hic_ch: GET_VALIDATED_CHANNEL(
                            "cram",
                            tolid_ver,
                            data.hic_data.hic_aligner,
                            data.assembly.defined_class,
                            data.assembly.project_id,
                            data.hic_data.hic_cram.collect()
                        )
            supplement_ch: tuple(
                [id: tolid_ver],
                data.assem_reads.supplementary_data,
            )
            genesets: (id == "FULL" || id == "JBROWSE" ? data.alignment.genesets.collect { geneset_path -> file(geneset_path, checkIfExists: true) } : [])
            synteny: (data.synteny ? data.synteny.collect { fasta -> file(fasta, checkIfExists: true) } : [])
            intron_size: (id == "FULL" ? data.intron.size : "")
            teloseq: data.telomere.teloseq
            busco_lineage: data.busco.lineage
            busco_lineages_path: file(data.busco.lineages_path, checkIfExists: true, type: 'dir')
        }
        .set { parsed }


    emit:
    ch_assembly_id    = parsed.tolid_version
    ch_reference      = parsed.reference
    ch_map_order      = parsed.map_order
    ch_assem_reads    = parsed.read_ch.filter { it } // filter []
    ch_kmer_prof_file = parsed.kmer_prof
    ch_hic_reads      = parsed.hic_ch
    ch_supp_reads     = parsed.supplement_ch
    ch_align_genesets = parsed.genesets.filter { it } // filter []
    ch_synteny_paths  = parsed.synteny.filter { it } // filter []
    ch_intron_size    = parsed.intron_size.filter { it } // filter ""
    ch_teloseq        = parsed.teloseq
    ch_lineageinfo    = parsed.busco_lineage
    ch_lineagespath   = parsed.busco_lineages_path
    versions          = ch_versions
}

def readYAML(yamlfile) {
    return new org.yaml.snakeyaml.Yaml().load(new FileReader(yamlfile.toString()))
}

def GET_VALIDATED_CHANNEL (data_type, tolid_ver, read_type, defined_class, project_id, data_files) {
    // Based on the the functions added in commit: 61f4ad9
    // Edited to be a function working on the raw yaml data
    // rather than channels as it was previously

    // Initialise defaults
    def fofn_files = []
    def direct_files = []
    def files_list = data_files instanceof List dePf? data_files : [data_files]

    // Process each file - separate FOFN from direct files
        files_list.each { file_path ->
            if (file_path.toString().contains(".fofn")) {
                def fofn_content = file(file_path).text.split('\n')
                    .collect { it.trim() }
                    .findAll { it } // Remove empty lines
                fofn_files.addAll(fofn_content)
            } else {
                direct_files.add(file_path)
            }
        }

        // Combine all files
        def all_files = direct_files + fofn_files

        // Validate files based on data type
        if (data_type == "cram") {
            def invalid_files = all_files.findAll {
                !it.toString().contains(".cram")
            }
            if (invalid_files.size() > 0) {
                error "One of the input hic files does not match cram format. Invalid files: ${invalid_files}"
            }
        } else if (data_type == "longread") {
            def invalid_files = all_files.findAll {
                !it.toString().contains(".fasta.gz") &&
                !it.toString().contains(".fa.gz") &&
                !it.toString().contains(".fn.gz")
            }
            if (invalid_files.size() > 0) {
                error "One of the input longread files does not match expected formats (fn.gz, fa.gz, fasta.gz). Invalid files: ${invalid_files}"
            }
        }

        // Create the resolved channel tuple
        def resolved_channel = tuple(
            [
                id:         tolid_ver,
                single_end: read_type != "illumina"
                aligner:    read_type && data_type == "cram" ? read_type : "NA",
                read_type:  read_type,
            ],
            all_files.collect { file(it, checkIfExists: true) }.unique()
        )

        return resolved_channel
    }
