---
title: 'TreeVal: A Nextflow pipeline for generating analyses to support manual curation of genome assemblies'
tags:
  - Python
  - Java
  - Nextflow
  - Groovy
  - Pipeline
  - Jbrowse2
  - HiGlass
  - Pretext
authors:
  - name: Ying Sims
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: William Eagles
    affiliation: 1
  - name: Damon-Lee Pointon
    affiliation: "1"
  - name: Shane McCarthy
    affiliation: "1"
affiliations:
 - name:
   index: 1
 - name:
   index: 2
date: 3 June 2024
bibliography: paper.bib

---

# Summary

The Tree of Life Project [@darwin2022sequence] aims to sequence the entirety of complex life in Britain and Northern Ireland, producing high-quality reference genome assemblies for an estimated 70,000 species in this region. Manual curation incorporates domain expertise to refine the output genome assemblies from the pipeline to a higher level of accuracy. We have developed a set of standardised analysis metrics that serve as guidelines for examining and evaluating genomes requiring curation. Additionally, we aim to automate the analysis process to achieve high-throughput efficiency. As a key step towards automation, we present the TreeVal data processing pipeline (v1.1.0), a bipartite project designed to generate and display a diverse array of data necessary for the manual curation of genomic assemblies into reference-quality assemblies. This is facilitated through tools such as JBrowse2 [@diesh2023jbrowse] and HiGlass [@kerpedjiev2018higlass]. The pipeline is implemented in Nextflow [@ewels2020nf] [@di2017nextflow], with dependencies distributed via conda [@anaconda], Docker [@merkel2014docker] and Singularity [@singularity]. The code is publicly available on GitHub at the following URL: https://github.com/sanger-tol/treeval/releases/tag/v1.1.0.


# Statement of need

The number of high-quality reference genome assemblies has increased dramatically in recent years with improvements to long-read sequencing technologies and the ramping up of Earth Biogenome Project (EBP) [@lewin2018earth] associated projects such as the Darwin Tree of Life Project [@darwin2022sequence], Aquatic Symbiosis Genomics Project [@mckenna2021aquatic] and Vertbrate Genomes Project [@rhie2021towards], among many others. Despite recent improvements to assembly and scaffolding algorithms taking advantage of these improved technologies, there is still a strong requirement for manual curation of these automated assemblies in order to convert them into chromosome-scale reference assemblies meeting or exceeding EBP standard metrics. To aid manual curation, it is desirable to prepare a standard set of analyses that can inform the curation process in making breaks and joins, separating haplotypes and assigning chromosomal units. Key components of the analysis include the generation of Hi-C contact maps, assessment of assembly read depth coverage, analysis of repeat element density, identification of sequence gaps, and annotation of telomere sequence. More in-depth genome structure analysis typically involves additional evaluations, such as self-comparative analysis to highlight regions of possible false duplication, transcriptome and proteome alignments, Benchmarking Universal Single-Copy Orthologs (BUSCO) analysis [@simao2015busco], and in silico digestion with Bionano restriction enzymes (optical maps). Historically, many of these analyses were coupled with Ensembl annotations, visualised through the gEVAL browser [@chow2016geval]. However, the gEVAL genome browser suffered from a lack of updates and was heavily reliant on a database-driven infrastructure, which posed significant maintenance challenges. As a result, we have transitioned to more flexible, flat-file-based viewers such as JBrowse2 [@diesh2023jbrowse], HiGlass [@kerpedjiev2018higlass] and PretextView [@edharry], which facilitate improved data management and provide domain experts with immediate visual feedback to evaluate genome assemblies. The production of these analyses has been implemented through the TreeVal pipeline, incorporating the aforementioned data processing methods and designed to be modular allowing analyses to be removed or new analyses added as methods evolve. The implementation of the TreeVal pipeline in Nextflow following nf-core [@ewels2020nf] [@di2017nextflow] guidelines will facilitate broader distribution of these analyses, enabling support for higher-throughput manual genome curation efforts across the wider EBP community. This pipeline has already supported the completion of over 1,000 reference genome assemblies from across the Tree of Life.

# Materials and Methods
## Features
This pipeline consists of sub-workflows which produce output files that are used to aid manual curation of the genome assembly, using the PretextView, HiGlass and JBrowse2 visualisation tools Diesh:2023. By using Nextflow DSL2 and following the nf-core guidelines, each singular process of the pipeline is encapsulated in a module and runs in an isolated environment (created by Conda [@anaconda], Docker [@merkel2014docker] or Singularity [@singularity]) allowing the workflow to be hardware agnostic. It also allows for the tuning of resources on a per process basis, if required, which is particularly useful when working on high performance compute clusters. In regards to these compute environments, nf-core also provides community created environment specific configuration files allowing the workflow to make use of the local systems queues and resources in a more dynamic manner.

## Workflow Overview

Each subworkflow comprising TreeVal (see Figure \autoref{fig:workflow}), is executed in parallel, when provided with sufficient resources. Treeval can be run in two modes, RAPID and FULL.

### Core Subworkflows (RAPID mode)

RAPID mode produces only the files required for the HiGlass and Pretext visualisation tools, so only uses a subset of the subworkflows. The **Hi-C Mapping** sub-workflow takes HiC read data and maps them against the input assembly. This has been implemented in what we have termed a super-module. Typically in the nf-core framework, only one to three tools are used per module to maximise portability. However, moving away from this has allowed for a significant saving in I/O and storage, which outweighs the cost of now being unable to finely tune the resources for each of the six tools used. These mapped segments, arising from the super-module, are then be merged so they can be ingested by PretextView (to generate a pretext graph, a modifiable Hi-C map) and produce an .mcool file for ingestion into HiGlass (along with a Juicer .hic file to display the HiC map in JBrowse2).

The remaining RAPID subworkflows produce files that can be ingested as tracks to annotate the HiC maps (Figure ???). The **Read Coverage** sub-workflow uses the PacBio long read data provided as .fasta.gz files, aligns them to the draft assembly using minimap2 [@li2018minimap] and outputs a .bigWig file that can be viewed as a coverage track. Additionally, the sub-workflow also produces punchlists for three different coverage depths: minimum, half, and maximum. This coverage information can be used in manual curation to detect regions of haplotypic duplication and in identifying sex chromosomes. The **Repeat Density** utilises WindowMasker [@Morgulis2006] and divides the genome into 10Kb bins, which act as sliding windows moving along the entire genome to profile the repeat intensity. Intersecting the bins with WindowMasker fragments and mapping the result produces a .bigWig file showing low-complexity repeat locations, which can be an indicator of potential centromeric regions. The **Gap Finder** subworkflow, employs seqtk cutn [@li2023seqtk] to cut any scaffold with a string of N’s. The .bed file is produced from calculating the lengths of these cut sequences. **Telo Finder** uses a user supplied telomeric sequence to identify putative telomeric regions in the draft assembly.

### Additional Subworkflows (FULL mode)

Running in FULL MODE produces tracks for the JBrowse2 tool, along with also running the RAPID subworkflows outlined above. The **In-silico Digest** sub-workflow is used to visualise the Bionano enzyme cutting sites for a genome .fasta file, generating JBrowse2 tracks based on three digestion enzymes (BSPQ1, BSSS1, and DLE1). For each of these, it identifies the recognition sequence of the labelling enzyme, mapping these to the reference genome .fasta producing an output in .bigBed format.

The two gene alignment sub-workflows utilises geneset data provided by the user and, ideally, closely related to the assembly being assesed. The **Nuclear Alignment** sub-workflow takes CDS, RNA and cDNA genesets, and aligns them against the input assembly with minimap2, outputting a .bigbed for display on JBrowse2. **Peptide Alignment** of he peptide genesets aligns using miniprot [@li2023miniprot] to produce a zipped .gff file. Both these methods, generate punchlists which list potential regions of interest (e.g., regions containing the same gene or a gene found over multiple regions). These can be opened in JBrowse using its spreadsheet viewer to aid the manual curation process.

The **BUSCO/Ancestral Gene** Analysis sub-workflow takes the draft assembly and extracts a list of BUSCO genes based on the v5 lineage ortholog gene database [@Manni2021] selected. Additionally, it provides an overlap BUSCO gene set based on a list of ancestral genes if available. Currently, these are only included for lepidoptera using genes identified corresponding to ancestral elements, known as Merian units, as determined by Wright [@Wright2024]. Both processes output their results for display as .bigBed files. Alternatively, the **Synteny** sub-workflow aligns provided syntenic genome files based on a selected clade. It then uses minimap2 to align these to the draft assembly using minimap2, producing a .paf for viewing.

Lastly, the self-comparative subworkflow (or **SelfComp**) employs four key approaches to generate the final alignment block result. First, it partitions the entire genome .fasta file into multiple 10kb units. This step is crucial to ensure a sufficient level of inter-alignment. The number of chunks is determined by the size of the genome, typically resulting in the fragmentation of a standard-sized genome (under 1Gb) into five portions. Second, the alignment process is executed between each pair of chunks using MUMMER [@marcais2018]. In our approach, the minimum match length is set at 400 base pairs, and it computes both forward and reverse complement matches. Additionally, it reports the query position of a reverse complement match relative to the forward strand of the query sequence. Third, the alignments obtained from MUMMER are associated with their original genomic locations. Finally, bedtools is used to concatenate alignments, allowing for 200Kb indels.

# Using the pipeline
## Configuration file
In order to execute the pipeline there must be a valid YAML file containing the locations of all files needed for the pipeline as well as values which act as modifiers to how the pipeline runs.

```yaml
assembly:
  assem_level: {scaffold|chromosome}
  assem_version: {Version number of the assembly}
  sample_id: {Name of the assembly, spaces should be replaced with underscores}
  latin_name: {Scientific Name}
  defined_class: {The user defined class of the input assembly}
  project_id: {Project ID of assembly} #Optional
reference_file: {Path to .f{a|n|asta}{.gz} formatted input genome }
map_order: length
assem_reads:
  read_type: hifi
  read_data: {Path to folder containing the longread data in .fasta.gz format}
hic_data:
  hic_cram: /lustre/scratch123/tol/resources/treeval/treeval-testdata/TreeValSmallData/Oscheius_DF5033/genomic_data/nxOscSpes1/hic-arima2/full/
  hic_aligner: {minimap2|bwamem2}
kmer_profile:
  kmer_length: {Default to 31}
  dir: {Path to pre-existing FKPROF files if they exist} #Optional
alignment:
  data_dir: {Path to the top level gene_alignment_data folder}
  geneset_id: {A csv delimited list of ScientificName.AssemblyName of the data to align}
self_comp:
  motif_len: 0
  mummer_chunk: 10
intron:
  size: "50k"
telomere:
  teloseq: {The expected telomeric sequence}
synteny:
  synteny_path: {Path to folder containing .fasta files used for syntenic alignments}
  synteny_genomes: {Specify the syntenic genome from the above path}
busco:
  lineages_path: {Path to the busco database e.g. /busco/v5}
  lineage: {The odb10 lineage to use}
```
\autoref{The TreeVal V1.1.0 input yaml format}

## Installation and Execution
The TreeVal pipeline contains three entry points specific to certain uses. The first does not require any explicit command and will run all subworkflows contained in the pipeline. We will refer to this as FULL, there is also RAPID and RAPID_TOL. The RAPID entry points execute a subset of the total subworkflows, focusing on generating files for visualisation in HiGlass (CITE) and PretextView (CITE). The difference between RAPID and RAPID_TOL is the pressense of a kmer plot generation subworkflow in the latter.

The pipeline can only be executed with docker or singulaity via the following command:
```
nextflow run sanger-tol/treeval -r 1.0.0 -profile {singularity|docker} --input {INPUT.yaml} --outdir {OUTDIR}
```
or
```
nextflow run sanger-tol/treeval -r 1.0.0 -profile {singularity|docker} --input {INPUT.yaml} --outdir {OUTDIR} -entry {RAPID|RAPID_TOL}
```
The pipeline can also be downloaded for offline use via the NF-core tools (CITE):
```
nf-core download sanger-tol/treeval -r 1.0.0
```

## Pre-processing of input data
More information on this topic can be found at the sanger pipelines website (https://pipelines.tol.sanger.ac.uk/treeval/1.1.0/usage) which is built from the `docs/usage.md` file of the TreeVal repository. It should be noted that there are some minimum requirements for data in this version of TreeVal. Pre-processed data should be stored in the following directory structure, missing data values:

```
treeval-resources
│
├─ gene_alignment_data/
│ └─ { defined_class }
│   ├─ csv_data
│   │ └─ { geneset_id }-data.csv # Generated by our scripts
│   └─ { geneset_id } # Here and below is generated by our scripts
│     └─ { geneset_id }
│       ├─ cdna
│       | └─ [Chunked fasta files]
│       ├─ rna
│       | └─ [Chunked fasta files]
│       ├─ cds
│       │ └─ [Chunked fasta files]
│       └─ peps
│         └─ [Chunked fasta files]
│
├─ synteny/       # Storage for your high quality genomes
| └─ { defined_class }
│   └─ [Syntenic Genomes .f{a|n|asta}{.gz}]
│
├─ treeval_yaml/  # Storage folder for you yaml files, it's useful to keep them
│
└─ treeval_stats/ # Storage for your TreeVal_run*.txt files
```
\autoref{The TreeVal V1.1.0 folder structure for pre-processed data}


1. The reference genome headers must be formatted to not include spaces or special characters, for example a header should be `>SCAFFOLD_1` rather than `>SCAFFOLD 1 @ CHROMOSOME 1`. Failure to comply will cause an error and halt the pipeline.

2. The geneset directory expects a specific data directory structure. An example would be as follows, values are taken from the TreeVal yaml file: `{alignment.data_dir}/{defined_class}/{csv_data}/{geneset_id}-data.csv`

This csv file contains the following information:
```
org,type,data_file
OsmiaBicornis.iOsmBic2,cds,/path/to/cds/OsmiaBicornis10004cds.MOD.fa
OsmiaBicornis.iOsmBic2,pep,/path/to/pep/OsmiaBicornis11021pep.MOD.fa
OsmiaBicornis.iOsmBic2,rna,/path/to/rna/OsmiaBicornis4001rna.MOD.fa
OsmiaBicornis.iOsmBic2,cdna,/path/to/cdna/OsmiaBicornis441cdna.MOD.fa
```

These geneset files are processed prior to pipeline execution with the Python3 (CITE) scripts found in the repository location `bin/treeval-dataprep/`. {NOTE: THE RUST ALTERNATIVE THAT REPLACES THIS IS PRETTY MUCH DONE - SHOULD WE JUST SWITCH NOW?}

3. Illumina HiC reads files should be in an unmapped CRAM format along side an index (`.cram.crai`) file.

4. Longread files should be in `fasta.gz` format. Most commonly, this will mean that data required reformatting from `fastq` format.

## Output

TreeVal produces a number of output files, many of which are generated for upload to JBrowse (CITE, these have been annotated with in brackets below). The `.pretext` are generated for use in PretextView, the contents of the `hic_files` folder (excluding `.pretext`) are used for visualisation in HiGlass.

```
{OUTDIR}
│
├─ treeval_upload/
| ├─ my.genome                # GENERATE_GENOME         (FULL)
│ ├─ coverage.bw              # READ_COVERAGE
│ ├─ coverage_log.bw          # READ_COVERAGE
| ├─ *_repeat_density.bw      # REPEAT_DENSITY
│ ├─ *gap.bed.gz              # GAP_FINDER
│ ├─ *gap.bed.gz.tbi          # GAP_FINDER
│ ├─ *telomere.bed.gz         # TELO_FINDER
│ ├─ *telomere.bed.gz.tbi     # TELO_FINDER
│ ├─ *_buscogene.bigbed       # BUSCO_ANALYSIS          (FULL)
│ ├─ *_ancestral.bigbed       # BUSCO_ANALYSIS          (FULL)
│ ├─ *.gff.gz                 # GENE_ALIGNMENT-PEPTIDE  (FULL)
│ ├─ *.gff.gz.tbi             # GENE_ALIGNMENT-PEPTIDE  (FULL)
│ ├─ *_cdna.bigBed            # GENE_ALIGNMENT-NUCLEAR  (FULL)
│ ├─ *_cds.bigBed             # GENE_ALIGNMENT-NUCLEAR  (FULL)
│ ├─ *_rna.bigBed             # GENE_ALIGNMENT-NUCLEAR  (FULL)
│ ├─ BSPQ1.bigBed             # INSILICO_DIGEST         (FULL)
│ ├─ BSSS1.bigBed             # INSILICO_DIGEST         (FULL)
│ ├─ DLE1.bigBed              # INSILICO_DIGEST         (FULL)
│ ├─ *_selfcomp.bigBed        # SELFCOMP                (FULL)
│ ├─ *.paf                    # SYNTENY                 (FULL)
│ ├─ *.ref.spectra-cn.ln.png  # KMER                    (FULL | RAPID_TOL)
│ ├─ *_{kmer_size}_.bw        # KMER_COVERAGE
│ └─ punchlists
│   ├─ halfcoverage.bigbed    # READ_COVERAGE
│   ├─ zerodepth.bigbed       # READ_COVERAGE
│   ├─ maxdepth.bigbed        # READ_COVERAGE
│   ├─ *_pep_punchlist.bed    # GENE_ALIGNMENT-PEPTIDE  (FULL)
│   ├─ *_cdna_punchlist.bed   # GENE_ALIGNMENT-NUCLEAR  (FULL)
│   ├─ *_cds_punchlist.bed    # GENE_ALIGNMENT-NUCLEAR  (FULL)
│   └─ *_rna_punchlist.bed    # GENE_ALIGNMENT-NUCLEAR  (FULL)
|
├─ hic_files
| ├─ *gap.bed                 # GAP_FINDER
│ ├─ *_repeat_density.bw      # REPEAT_DENSITY
│ ├─ *.mcool                  # HIC_MAPPING
│ ├─ *_pretext_normal.pretext # HIC_MAPPING + PRETEXT_INGESTION
│ ├─ *_pretext_hr.pretext     # HIC_MAPPING + PRETEXT_INGESTION
│ └─ *telomere.bed            # TELO_FINDER
|
└─ pipeline_info
  ├─ TreeVal_Runs*.txt        # TreeValProject.Summary (Custom Groovy Function)
  ├─ execution*{html|txt}     # STANDARD OUTPUT
  ├─ pipeline*{html|txt}      # STANDARD OUTPUT
  └─ software_versions.yml    # STANDARD OUTPUT

```
\autoref{The TreeVal V1.1.0 output directory structure}

# Conclusions and Discussions
The TreeVal pipeline provides a comprehensive set of functionalities to facilitate genome assembly manual curation in a single execution. Users can avoid the hassle of downloading dependencies and installing relevant packages and environments, thanks to the user-friendly Nextflow framework [@ewels2020nf] [@di2017nextflow]. We are currently investigating ways to scale up the analysis pipeline, emphasizing flexibility and dynamism. In this effort, we plan to deliver a more versatile configuration profile that accommodates genomes of various sizes and diverse genomic characteristics. Additionally, we are implementing functions that accept optional arguments, enabling the selection of specific subworkflows for execution. Furthermore, we are exploring additional analyses such as linkage analysis and further grouping of ancestral elements. Despite these ongoing enhancements, we believe this pipeline is a robust platform that supports the goals of the Darwin Tree of Life (DTOL) project [@darwin2022sequence] and the broader Earth BioGenome Project (EBP) [@lewin2018earth].

# Figures
![An Overview of the TreeVal Workflow, showing component subworkflows and their outputs. Also shown are the suggested visualisation tools the outputs can be viewed in.](./workflow_overview.png)

# Acknowledgements
To advance this project, we acknowledge the contributions of the Genome Reference Informatics Team from The Wellcome Sanger Institute Tree of Life programme for their valuable advice on analysis methods and validation of pipeline outcomes. We extend our gratitude to all members of the Tree of Life Assembly Team and The Tree of Life IT Team for providing essential infrastructure support, data storage, and management services. Additionally, our appreciation goes out to the Nextflow and NF-Core communities for their collaboration and support.

# References