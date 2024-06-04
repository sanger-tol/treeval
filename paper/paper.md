---
title: 'TreeVal: A Nextflow pipeline for generating analyses to support manual curation of genome assemblies'
tags:
  - Python
  - Java
  - Nextflow
authors:
  - name: Ying Sims
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Will Eagle
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

The Tree of Life Project [@darwin2022sequence] aims to sequence the entirety of complex life in Britain and Northern Ireland, producing high-quality reference genome assemblies for an estimated 70,000 species in this region. Manual curation incorporates domain expertise to refine the output genome assemblies from the pipeline to a higher level of accuracy. We have developed a set of standardised analysis metrics that serve as guidelines for examining and evaluating genomes requiring curation. Additionally, we aim to automate the analysis process to achieve high-throughput efficiency. As a key step towards automation, we present the TreeVal data processing pipeline (v1.1.0), a bipartite project designed to generate and display a diverse array of data necessary for the manual curation of genomic assemblies into reference-quality assemblies. This is facilitated through tools such as JBrowse2 [@diesh2023jbrowse] and HiGlass [@kerpedjiev2018higlass]. The pipeline is implemented in Nextflow [@ewels2020nf] [@di2017nextflow], with dependencies distributed via conda, Docker and Singularity. The code is publicly available on GitHub at the following URL: https://github.com/sanger-tol/treeval/releases/tag/v1.1.0.


# Statement of need

The number of high-quality reference genome assemblies has increased dramatically in recent years with improvements to long-read sequencing technologies and the ramping up of Earth Biogenome Project (EBP) [@lewin2018earth] associated projects such as the Darwin Tree of Life Project [@darwin2022sequence], Aquatic Symbiosis Genomics Project [@mckenna2021aquatic] and Vertbrate Genomes Project [@rhie2021towards], among many others. Despite recent improvements to assembly and scaffolding algorithms taking advantage of these improved technologies, there is still a strong requirement for manual curation of these automated assemblies in order to convert them into chromosome-scale reference assemblies meeting or exceeding EBP standard metrics. To aid manual curation, it is desirable to prepare a standard set of analyses that can inform the curation process in making breaks and joins, separating haplotypes and assigning chromosomal units. Key components of the analysis include the generation of Hi-C contact maps, assessment of assembly read depth coverage, analysis of repeat element density, identification of sequence gaps, and annotation of telomere sequence. More in-depth genome structure analysis typically involves additional evaluations, such as self-comparative analysis to highlight regions of possible false duplication, transcriptome and proteome alignments, Benchmarking Universal Single-Copy Orthologs (BUSCO) analysis [@simao2015busco], and in silico digestion with Bionano restriction enzymes (optical maps). Historically, many of these analyses were coupled with Ensembl annotations, visualised through the gEVAL browser [@chow2016geval]. However, the gEVAL genome browser suffered from a lack of updates and was heavily reliant on a database-driven infrastructure, which posed significant maintenance challenges. As a result, we have transitioned to more flexible, flat-file-based viewers such as JBrowse2 [@diesh2023jbrowse], HiGlass [@kerpedjiev2018higlass] and PretextView [@edharry], which facilitate improved data management and provide domain experts with immediate visual feedback to evaluate genome assemblies. The production of these analyses has been implemented through the TreeVal pipeline, incorporating the aforementioned data processing methods and designed to be modular allowing analyses to be removed or new analyses added as methods evolve. The implementation of the TreeVal pipeline in Nextflow following nf-core [@ewels2020nf] [@di2017nextflow] guidelines will facilitate broader distribution of these analyses, enabling support for higher-throughput manual genome curation efforts across the wider EBP community. This pipeline has already supported the completion of over 1,000 reference genome assemblies from across the Tree of Life.

# Materials and Methods
## Features


## Workflow
Figure 1 \autoref{workflow_overview} depicts the workflow:



## Configuration file


## Installation and Execution


## Output and Error Handling


# Conclusions and Discussions
The TreeVal pipeline provides a comprehensive set of functionalities to facilitate genome assembly manual curation in a single execution. Users can avoid the hassle of downloading dependencies and installing relevant packages and environments, thanks to the user-friendly Nextflow framework [@ewels2020nf] [@di2017nextflow]. We are currently investigating ways to scale up the analysis pipeline, emphasizing flexibility and dynamism. In this effort, we plan to deliver a more versatile configuration profile that accommodates genomes of various sizes and diverse genomic characteristics. Additionally, we are implementing functions that accept optional arguments, enabling the selection of specific subworkflows for execution. Furthermore, we are exploring additional analyses such as linkage analysis and further grouping of ancestral elements. Despite these ongoing enhancements, we believe this pipeline is a robust platform that supports the goals of the Darwin Tree of Life (DTOL) project [@darwin2022sequence] and the broader Earth BioGenome Project (EBP) [@lewin2018earth].


# Figures
![An Overview of the TreeVal Workflow, showing component subworkflows and their outputs. Also shown are the suggested visualisation tools the outputs can be viewed in.](./workflow_overview.png)

# Acknowledgements


# References