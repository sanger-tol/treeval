# sanger-tol/treeval: Output

# Introduction

This document describes the output produced by the TreeVal pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

# Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following workflows:

- [generate-genome](#generate-genome) - Builds a genome description file of the reference genome.
- [read-coverage](#read-coverage) - Produces read coverage based on hifi/clr/ont/illumina fasta file/s.
- [gap-finder](#gap-finder) - Identifies gaps in the input genome.
- [repeat-density](#repeat-density) - Reports the intensity of regional repeats within an input assembly.
- [hic-mapping](#hic-mapping) - Aligns illumina HiC short reads to the input genome, generates mapping file in three format for visualisation: `.pretext`, `.hic` and `.mcool`.
- [telo-finder](#telo-finder) - Identifies regions of a user given telomeric sequence.
- [gene-alignment](#gene-alignment) - Aligns the peptide and nuclear data from assemblies of related species to the input genome.
- [insilico-digest](#insilico-digest) - Generates a map of enzymatic digests using 3 Bionano enzymes.
- [selfcomp](#selfcomp) - Identifies regions of self-complementary sequence.
- [synteny](#synteny) - Generates syntenic alignments between the input and other high quality genomes.
- [busco-analysis](#busco-analysis) - Uses BUSCO to identify ancestral elements. Also use to identify ancestral Lepidopteran genes (merian units).
- [kmer](#kmer) - Counts k-mer and generates a copy number spectra plot.
- [kmer coverage](#kmer-coverage) - Counts k-mer (or uses existing k-mer profile) and produces k-mer coverage.
- [pretext-ingestion](#pretext-ingestion) - Ingests accessory files into the pretext file.

- [pipeline-information](#pipeline-information) - Report metrics generated during the workflow execution

## generate-genome

This workflow generates a .genome file which describes the base pair length of each scaffold in the reference genome. This file is then recycled into the workflow to be used by a number of other subworkflows.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `my.genome`: Genome description file of the reference genome.

</details>

![Generate genome workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_generate_genome.png)

## read-coverage

Read Coverage uses genome sequence reads (HiFi, CLR, ONT or Illumina) reads to generate a coverage bigWig as well as a trio of depth.bigbed files.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `coverage.bw`: Coverage of aligned reads across the reference genome in bigwig format.
  - `coverage_log.bw`: A log corrected coverage file which aims to smooth out the above track.
- `treeval_upload/punchlists/`
  - `maxdepth.bigbed`: Max read depth punchlist in bigBed format.
  - `zerodepth.bigbed`: Zero read depth punchlist in bigBed format.
  - `halfcoverage.bigbed`: Half read depth punchlist in bigBed format.

</details>

![Read Coverage workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_read_coverage.png)

## gap-finder

The gap-finder subworkflow generates a bed file containing the genomic locations of the gaps in the sequence. This file is injected into the hic_maps at a later stage. The output bed file is then BGzipped and indexed for display on JBrowse.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`: A bgzipped file containing gap locations
  - `*.bed.gz.tbi`: A tabix index file for the above file.
- `hic_files/`
  - `*.bed`: The raw bed file needed for ingestion into Pretext

</details>

![Gap Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_gap_finder.png)

## repeat-density

This uses [WindowMasker](https://github.com/goeckslab/WindowMasker) to mark potential repeats on the genome. The genome is chunked into 10kb bins which move along the entire genome as sliding windows in order to profile the repeat intensity. These fragments are then mapped back to the original assembly for visualisation purposes.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_repeat_density.bw`: Intersected read windows aligned to the reference genome in bigwig format.

</details>

![Repeat Density workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_repeat_density.png)

## hic-mapping

The hic-mapping subworkflow takes a set of HiC read files in .cram format as input and derives HiC mapping outputs in .pretext, .hic, and .mcool formats. These outputs are used for visualization on [PretextView](https://github.com/wtsi-hpag/PretextView), [Juicebox](https://github.com/aidenlab/Juicebox), and [HiGlass](https://github.com/higlass/higlass) respectively.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_pretext_hr.pretext`: High resolution pretext map.
  - `*_pretext_normal.pretext`: Standard resolution pretext map.
  - `*.mcool`: HiC map required for HiGlass

</details>

![Hic Mapping workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_hic_mapping.png)

## telo-finder

The telo-finder subworkflow uses a supplied (by the .yaml) telomeric sequence to identify putative telomeric regions in the input genome. The BGZipped and indexed file is used in JBrowse and as supplementary data for HiGlass and PreText.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`: A bgzipped file containing telomere sequence locations
  - `*.bed.gz.tbi`: A tabix index file for the above file.
- `hic_files/`
  - `*.bed`: The raw .bed file needed for ingestion into Pretext

</details>

![Telomere Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_telo_finder.png)

## busco-analysis

The BUSCO annotation subworkflow takes an assembly genome as input and extracts a list of [BUSCO](https://gitlab.com/ezlab/busco) genes based on the BUSCO results obtained from BUSCO. Additionally, it provides an overlap BUSCO gene set based on a list of lepidoptera ancestral genes (Wright et al., 2023), which has been investigated by Charlotte Wright from Mark Blaxter's lab at the Sanger Institute.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_buscogene.bigbed`: BigBed file for BUSCO genes track.
  - `*_ancestral.bigbed`: BigBed file for ancestral elements track.

</details>

![Busco analysis workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_busco_analysis.png)

## gene-alignment

The gene alignment subworkflows load genesets (cdna, cds, rna, pep) data from a csv list of genomes, in the input .yaml, and aligns these to the reference genome. It contains two subworkflows, one of which handles peptide data and the other of which handles RNA, CDS and complementary DNA data. These produce files that can be displayed by JBrowse as tracks.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.gff.gz`: Zipped .gff for each species with peptide data.
  - `*.gff.gz.tbi`: TBI index file of each zipped .gff.
  - `*_cdna.bigBed`: BigBed file for each species with complementary DNA data.
  - `*_cds.bigBed`: BigBed file for each species with nuclear DNA data.
  - `*_rna.bigBed`: BigBed file for each species with nRNA data.
- `treeval_upload/punchlists/`
  - `*_pep_punchlist.bed`: Punchlist for peptide track.
  - `*_cdna_punchlist.bed`: Punchlist for cdna track.
  - `*_cds_punchlist.bed`: Punchlist for cds track.
  - `*_rna_punchlist.bed`: Punchlist for rna track.

</details>

![Gene alignment workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_gene_alignment.png)

## insilico-digest

The insilico-digest workflow is used to visualize the Bionano enzyme cutting sites for a genomic FASTA file. This procedure generates data tracks based on three digestion enzymes (by default): BSPQ1, BSSS1, and DLE1.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `{BSPQ1|BSSS1|DLE1}.bigBed`: Bionano insilico digest cut sites track in the bigBed format for each of the set digestion enzymes.

</details>

![Insilico digest workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_insilico_digest.png)

## selfcomp

The selfcomp subworkflow is a comparative genomics analysis algorithm originally performed by the Ensembl projects database, and reverse engineered in Python3 by @yumisims. It involves comparing the genes and genomic sequences within a single species. The goal of the analysis is to identify haplotypic duplications in a particular genome assembly.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_selfcomp.bigBed`: BigBed file containing selfcomp track data.

</details>

![Selfcomp workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_self_comp.png)

## synteny

This subworkflow searches along a predetermined path for syntenic genome files based on clade and then aligns with [MINIMAP2_ALIGN](https://nf-co.re/modules/minimap2_align) to the reference genome, emitting an aligned .paf file for each.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.paf`: .paf file for each syntenic genomic aligned to reference.

</details>

![Synteny workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_synteny.png)

## kmer

This subworkflow performs a k-mer count using [FASTK_FASTK](https://nf-co.re/modules/fastk_fastk) then passes the results to [MERQURYFK_MERQURYFK](https://nf-co.re/modules/merquryfk_merquryfk) to plot a copy-number k-mer spectra.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*.ref.spectra-cn.ln.png`: .png file of copy number k-mer spectra.

</details>

![Kmer Workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_kmer.png)

## kmer coverage

This subworkflow performs a k-mer count using [FASTK_FASTK](https://nf-co.re/modules/fastk_fastk) (or uses an already existing k-mer profile) then passes the results to FKUTILS_FKPROF to produces k-mer coverage track.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_{kmer_size}_.bw`: .png file of copy number k-mer spectra.

</details>

![Kmer coverage Workflow](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_kmer_coverage.png)

## Full Workflow diagram

The full pipeline diagram is very large, with the pipeline consisting of over 100 processes.
![The Pipeline](https://raw.githubusercontent.com/sanger-tol/treeval/main/docs/images/v1-1-0/treeval_1_1_0_full_diagram.png)

## pipeline-information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `treeval_info/`
  - `TreeVal_Runs*.txt`: A Report generated by TreeValProject.
  - `*txt`: Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>
