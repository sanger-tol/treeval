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
            read_ch: (data.assem_reads.read_type in ['hifi', 'clr', 'ont', 'illumina']
                ? tuple(
                    [
                        id: tolid_ver,
                        single_end: data.assem_reads.read_type != "illumina",
                        read_type: data.assem_reads.read_type,
                    ],
                    data.assem_reads.read_data.collect { fq -> file(fq, checkIfExists: true, type: 'dir') },
                )
                : [])
            kmer_prof: tuple(
                [
                    id: tolid_ver,
                    kmer: kmer_len,
                ],
                file("${data.kmer_profile.dir}/k${kmer_len}/${data.assembly.sample_id}.k${kmer_len}.ktab"),
            )
            hic_ch: tuple(
                [
                    id: tolid_ver,
                    aligner: data.hic_data.hic_aligner,
                ],
                data.hic_data.hic_cram.collect { cram -> file(cram, checkIfExists: true, type: 'dir')},
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
