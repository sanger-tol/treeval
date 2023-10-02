# sanger-tol/treeval: Output

# Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

# Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following workflows:

- [generate-genome](#generate-genome) - Builds genome description file of the reference genome.
- [longread-coverage](#longread-coverage) - Produces read coverage based on pacbio long read fasta file.
- [gap-finder](#gap-finder) - Identifies contig gaps in the input genome.
- [repeat-density](#repeat-density) - Reports the intensity of regional repeats within an input assembly.
- [hic-mapping](#hic-mapping) - Aligns illumina HiC short reads to the input genome, generates mapping file in three format for visualisation: .pretext, .hic and .mcool
- [telo-finder](#telo-finder) - Identifies regions of a user given telomeric sequence.
- [gene-alignment](#gene-alignment) - Aligns the peptide and nuclear data from assemblies of related species to the input genome.
- [insilico-digest](#insilico-digest) - Generates a map of enzymatic digests using 3 Bionano enzymes.
- [selfcomp](#selfcomp) - Identifies regions of self-complementary sequence.
- [synteny](#synteny) - Generates syntenic alignments between other high quality genomes.
- [busco-analysis](#busco-analysis) - Uses BUSCO to identify ancestral elements. Also use to identify ancestral Lepidopteran genes (merian units).

- [pipeline-information](#pipeline-information) - Report metrics generated during the workflow execution

## generate-genome

This workflow generates a .genome file which describes the base pair length of each scaffold in the reference genome. This file is then recycled into the workflow to be used by a number of other subworkflows.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `my.genome`: Genome description file of the reference genome.

</details>

![Generate genome workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_generate_genome.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## longread-coverage

Longread Coverage uses Pacbio HiC reads to generage a coverage bigWig as well as a trio of depth.bigbed files.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `coverage.bw`: Coverage of aligned reads across the reference genome in bigwig format.
- `treeval_upload/punchlists/`
  - `maxdepth.bigbed`: Max read depth punchlist in bigBed format.
  - `zerodepth.bigbed`: Zero read depth punchlist in bigBed format.
  - `halfcoverage.bigbed`: Half read depth punchlist in bigBed format.

</details>

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

![Gap Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_gap_finder.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## repeat-density

This uses [WindowMasker](https://github.com/goeckslab/WindowMasker) to mark potential repeats on the genome. The genome is chunked into 10kb bins which move along the entire genome as sliding windows in order to profile the repeat intensity. These fragments are then mapped back to the original assembly for visualization purposes.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_repeat_density.bw`: Intersected read windows aligned to the reference genome in bigwig format.

</details>

![Repeat Density workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_repeat_density.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## hic-mapping

The hic-mapping subworkflow takes a set of HiC read files in .cram format as input and derives HiC mapping outputs in .pretext, .hic, and .mcool formats. These outputs are used for visualization on [PretextView](https://github.com/wtsi-hpag/PretextView), [Juicebox](https://github.com/aidenlab/Juicebox), and [Higlass](https://github.com/higlass/higlass) respectively.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_pretext_hr.pretext`: High resolution pretext map.
  - `*_pretext_normal.pretext`: Low resolution pretext map.
  - `*.mcool`: HiC map required for HiGlass

</details>

![Hic Mapping workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_hic_mapping.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

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

![Telomere Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_telo_finder.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## busco-analysis

The BUSCO annotation subworkflow takes an assembly genome as input and extracts a list of [BUSCO](https://gitlab.com/ezlab/busco) genes based on the BUSCO results obtained from BUSCO. Additionally, it provides an overlap BUSCO gene set based on a list of lepidoptera ancestral genes (Wright et al., 2023), which has been investigated by Charlotte Wright from Mark Blaxter's lab at the Sanger Institute.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_buscogene.bigbed`: BigBed file for BUSCO genes track.
  - `*_ancestral.bigbed`: BigBed file for ancestral elements track.

</details>

![Busco analysis workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_busco_analysis.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## gene-alignment

The gene alignment subworkflows load genesets (cdna, cds, rna, pep) data from a given list of genomes detailed, in the input .yaml, and aligns these to the reference genome. It contains two subworkflows, one of which handles peptide data and the other of which handles RNA, nuclear and complementary DNA data. These produce files that can be displayed by JBrowse as tracks.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.gff.gz`: Zipped .gff for each species with peptide data.
  - `*.gff.gz.tbi`: TBI index file of each zipped .gff.
  - `*_cdna.bigBed`: BigBed file for each species with complementary DNA data.
  - `*_cds.bigBed`: BigBed file for each species with nuclear DNA data.
  - `*_rna.bigBed`: BigBed file for each species with nRNAdata.
- `treeval_upload/punchlists/`
  - `*_pep_punchlist.bed`: Punchlist for peptide track.
  - `*_cdna_punchlist.bed`: Punchlist for cdna track.
  - `*_cds_punchlist.bed`: Punchlist for cds track.
  - `*_rna_punchlist.bed`: Punchlist for rna track.

</details>

![Gene alignment workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_gene_alignment.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## insilico-digest

The insilico-digest workflow is used to visualize the Bionano enzyme cutting sites for a genomic FASTA file. This procedure generates data tracks based on three digestion enzymes: BSPQ1, BSSS1, and DLE1.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `{BSPQI|BSSSI|DLE1}.bigBed`: Bionano insilico digest cut sites track in the bigBed format for each of the set digestion enzymes.

</details>

![Insilico digest workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_insilico_digest.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## selfcomp

The selfcomp subworkflow is a comparative genomics analysis originally performed by the Ensembl project. It involves comparing the genes and genomic sequences within a single species. The goal of the analysis is mainly to identify haplotypic duplications in a particular genome assembly.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_selfcomp.bigBed`: BigBed file containing selfcomp track data.

</details>

![Selfcomp workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_selfcomp.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## synteny

This worflows searches along predetermined path for syntenic genome files based on clade and then aligns with [MINIMAP2_ALIGN](https://nf-co.re/modules/minimap2_align) each to the reference genome, emitting an aligned .paf file for each.

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.paf`: .paf file for each syntenic genomic aligned to reference.

</details>

![Synteny workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_synteny.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## kmer

This worflows performs a k-mer count using [FASTK_FASTK](https://nf-co.re/modules/fastk_fastk) then passes the results to [MERQURYFK_MERQURYFK](https://nf-co.re/modules/merquryfk_merquryfk) to plot a copy-number k-mer spectra.

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*.ref.spectra-cn.ln.png`: .png file of copy number k-mer spectra.

</details>

![Synteny workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_synteny.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

## pipeline-information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `treeval_info/`
  - Report generated by TreeValProject.
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>
