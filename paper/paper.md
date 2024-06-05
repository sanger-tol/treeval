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
This pipeline consists of sub-workflows that produce output files used to aid manual curation of the genome assembly, utilizing the PretextView, HiGlass, and JBrowse2 visualization tools [@diesh2023jbrowse]. By using Nextflow DSL2 and following the nf-core guidelines, each process of the pipeline is encapsulated in a module which runs in an isolated environment (In the case of TreeVal V1.1.0 this environment is generated by Docker [@merkel2014docker] or Singularity [@singularity]) allowing the workflow to be hardware agnostic. It also allows for the tuning of resources on a per process basis, if required, which is particularly useful when working on high performance compute clusters with diverse data. In regards to these compute environments, nf-core also provides community created environment specific configuration files allowing the workflow to make use of the local systems queues and resources in a more dynamic manner.

## Workflow Overview

Each subworkflow comprising TreeVal (see Figure \autoref{fig:workflow}), is executed in parallel, when provided with sufficient resources. Treeval can be run in two modes, RAPID and FULL.

### Core Subworkflows (RAPID mode)

RAPID mode produces only the files required for the HiGlass and Pretext visualisation tools, so only uses a subset of the subworkflows. The **Hi-C Mapping** sub-workflow takes HiC read data and maps them against the input assembly. This has been implemented in what we have termed a super-module. Typically in the nf-core framework, only one to three tools are used per module to maximise portability. However, moving away from this has allowed for a significant saving in I/O and storage, which outweighs the cost of now being unable to finely tune the resources for each of the six tools used. These mapped segments, arising from the super-module, are then be merged so they can be ingested by PretextView (to generate a pretext graph, a modifiable Hi-C map) and produce an .mcool file for ingestion into HiGlass (along with a Juicer .hic file to display the HiC map in JBrowse2).

The remaining RAPID subworkflows produce files that can be ingested as tracks to annotate the HiC maps (Figure ???). The **Read Coverage** sub-workflow uses the PacBio long read data provided as .fasta.gz files, aligns them to the draft assembly using minimap2 [@li2018minimap] and outputs a .bigWig file that can be viewed as a coverage track. Additionally, the sub-workflow also produces punchlists for three different coverage depths: minimum, half, and maximum. This coverage information can be used in manual curation to detect regions of haplotypic duplication and in identifying sex chromosomes. The **Repeat Density** utilises WindowMasker [@Morgulis2006] and divides the genome into 10Kb bins, which act as sliding windows moving along the entire genome to profile the repeat intensity. Intersecting the bins with WindowMasker fragments and mapping the result produces a .bigWig file showing low-complexity repeat locations, which can be an indicator of potential centromeric regions. The **Gap Finder** subworkflow, employs seqtk cutn [@li2023seqtk] to cut any scaffold with a string of N’s. The .bed file is produced from calculating the lengths of these cut sequences. **Telo Finder** uses a user supplied telomeric sequence to identify putative telomeric regions in the draft assemble

### Additional Subworkflows (FULL mode)

Running in FULL mode produces tracks for the JBrowse2 genome browser, along with also running the RAPID subworkflows outlined above. The **In-silico Digest** sub-workflow is used to visualise the Bionano enzyme cutting sites for a genomic FASTA file, generating JBrowse2 tracks based on three digestion enzymes (BSPQ1, BSSS1, and DLE1). For each of these, it identifies the recognition sequence of the labelling enzyme, mapping these to the reference genome .fasta producing an output in .bigBed format.

The two gene alignment sub-workflows utilises geneset data provided by the user and, ideally, closely related to the assembly being assesed. The **Nuclear Alignment** sub-workflow takes CDS, RNA and cDNA genesets, and aligns them against the input assembly with minimap2, outputting a .bigbed for display on JBrowse2. **Peptide Alignment** of he peptide genesets aligns using miniprot [@li2023miniprot] to produce a zipped .gff file. Both these methods, generate punchlists which list potential regions of interest (e.g., regions containing the same gene or a gene found over multiple regions). These can be opened in JBrowse using its spreadsheet viewer to aid the manual curation process.

The **BUSCO/Ancestral Gene Analysis** sub-workflow takes the draft assembly and extracts a list of BUSCO genes based on the v5 lineage ortholog gene database [@Manni2021] selected. Additionally, it provides an overlap BUSCO gene set based on a list of ancestral genes if available. Currently, these are only included for lepidoptera using genes identified corresponding to ancestral elements, known as Merian units, as determined by Wright @Wright2024. Both processes output their results for display as .bigBed files. Alternatively, the **Synteny** sub-workflow aligns provided syntenic genome files based on a selected clade. It then uses minimap2 to align these to the draft assembly using minimap2, producing a .paf for viewing.

Lastly, the self-complementary subworkflow (or **SelfComp**) employs four key approaches to generate the final alignment block result. First, it partitions the entire input FASTA file into multiple 10kb units. This step is crucial to ensure a sufficient level of inter-alignment. The number of chunks is determined by the size of the genome, typically resulting in the fragmentation of a standard-sized genome (under 1Gb) into five portions. Second, the alignment process is executed between each pair of chunks using MUMMER [@marcais2018]. In our approach, the minimum match length is set at 400 base pairs, and it computes both forward and reverse complement matches. Additionally, it reports the query position of a reverse complement match relative to the forward strand of the query sequence. Third, the alignments obtained from MUMMER are associated with their original genomic locations. Finally, bedtools is used to concatenate alignments, allowing for 200Kb indels.

# Using the pipeline
## Configuration and Execution
In order to execute the pipeline there must be a valid yaml file containing the locations of all files needed for the pipeline as well as values which act as modifiers to how the pipeline runs.

CHEAT SHEET - containing

## Output

TreeVal produces a number of output files, many of which are generated for upload to JBrowse2 (these have been annotated with in brackets below). The `.pretext` are generated for use in PretextView, the contents of the `hic_files` folder (excluding `.pretext`) are used for visualisation in HiGlass.

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
│ ├─ *.ref.spectra-cn.ln.png  # KMER                    (FULL)
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
The TreeVal pipeline provides a comprehensive set of functionalities to facilitate genome assembly and manual curation in a single execution. Users can avoid the hassle of downloading dependencies and installing relevant packages and environments, thanks to the user-friendly Nextflow framework [@ewels2020nf] [@di2017nextflow]. We are currently investigating ways to scale up the analysis pipeline, emphasizing flexibility and dynamism. In this effort, we plan to deliver a more versatile configuration profile that accommodates genomes of various sizes and diverse genomic characteristics. Additionally, we are implementing functions that accept optional arguments, enabling the selection of specific subworkflows for execution. Furthermore, we are exploring additional analyses such as linkage analysis and further grouping of ancestral elements. Despite these ongoing enhancements, we believe this pipeline is a robust platform that supports the goals of the Darwin Tree of Life (DTOL) project [@darwin2022sequence] and the broader Earth BioGenome Project (EBP) [@lewin2018earth].

# Figures
![An Overview of the TreeVal Workflow, showing component subworkflows and their outputs. Also shown are the suggested visualisation tools the outputs can be viewed in.](./workflow_overview.png)

# Acknowledgements
To advance this project, we acknowledge the contributions of the Genome Reference Informatics Team from The Wellcome Sanger Institute Tree of Life programme for their valuable advice on analysis methods and validation of pipeline outcomes. We extend our gratitude to all members of the Tree of Life Assembly Team and The Tree of Life IT Team for providing essential infrastructure support, data storage, and management services. Additionally, our appreciation goes out to the Nextflow and NF-Core communities for their collaboration and support.

# References