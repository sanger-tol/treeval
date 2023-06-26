# sanger-tol/treeval: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following workflows:

- [YAML_INPUT](#YAML_INPUT) - Reads the input yaml and generates parameters used by other workflows.
- [GENERATE_GENOME](#GENERATE_GENOME) - Builds genome description file of the reference genome.
- [LONGREAD_COVERAGE](#LONGREAD_COVERAGE) - Produces read coverage based on pacbio long read fasta file.
- [GAP_FINDER](#GAP_FINDER) - Identifies contig gaps in the input genome.
- [REPEAT_DENSITY](#REPEAT_DENSITY) - Reports the intensity of regional repeats within an input assembly.
- [HIC_MAPPING](#HIC_MAPPING) - Aligns illumina HiC short reads to the input genome, generates mapping file in three format for visualisation: .pretext, .hic and .mcool
- [TELO_FINDER](#TELO_FINDER) - Identifies regions of a user given telomeric sequence.
- [GENE_ALIGNMENT](#GENE_ALIGNMENT) - Aligns the peptide and nuclear data from assemblies of related species to the input genome.
- [INSILICO_DIGEST](#INSILICO_DIGEST) - Generates a map of enzymatic digests using 3 Bionano enzymes.
- [SELFCOMP](#SELFCOMP) - Identifies regions of self-complementary sequence.
- [SYNTENY](#SYNTENY) - Generates syntenic alignments between other high quality genomes.
- [BUSCO_ANALYSIS](#BUSCO_ANALYSIS) - Uses BUSCO to identify ancestral elements. Also use to identify ancestral Lepidopteran genes (merian units).

- [Pipeline information](#Pipeline-information) - Report metrics generated during the workflow execution

### YAML_INPUT

This subworkflow reads the input .yaml via the use of the built-in snakeyaml.Yaml component, which converts the yaml into a nested list. Via some simple channel manipulation, each item in this nexted list is converted into a parameter for use in each of the other subworkflows.

### GENERATE_GENOME

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `my.genome`: Genome description file of the reference genome.

</details>

This workflow generates a .genome file which describes the base pair length of each scaffold in the reference genome. This is performed by [SAMTOOLS_FAIDX](https://nf-co.re/modules/samtools_faidx) to generate a .fai file. This index file is trimmed using local module [GENERATE_GENOME_FILE](../modules/local/generate_genome_file.nf) to output a .genome file. This file is then recycled into the workflow to be used by a number of other subworkflows.

![Generate genome workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_generate_genome.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### LONGREAD_COVERAGE

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `coverage.bw`: Coverage of aligned reads across the reference genome in bigwig format.
- `treeval_upload/punchlists/`
  - `maxdepth.bigbed`: Max read depth punchlist in bigBed format.
  - `zerodepth.bigbed`: Zero read depth punchlist in bigBed format.
  - `halfcoverage.bigbed`: Half read depth punchlist in bigBed format.

</details>

Longread Coverage uses Pacbio HiC reads to generage a coverage bigWig as well as a trio of depth.bigbed files

### GAP_FINDER

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`: A bgzipped file containing gap locations
  - `*.bed.gz.tbi`: A tabix index file for the above file.
- `hic_files/`
  - `*.bed`: The raw bed file needed for ingestion into Pretext

</details>

The GAP_FINDER subworkflow generates a bed file containing the genomic locations of the gaps in the sequence. This is performed by the use of [SEQTK_CUTN]() which cuts the input genome at sites of N (gaps). [GAP_LENGTH]() then calculates the lengths of gaps generates in the previous step, this file is injected into the hic_maps at a later stage. SEQTK's output bed file is then BGzipped and indexed by [TABIX_BGZIPTABIX](https://nf-co.re/modules/tabix_bgziptabix).

![Gap Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_gap_finder.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### REPEAT_DENSITY

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_repeat_density.bw`: Intersected read windows aligned to the reference genome in bigwig format.

</details>
This uses [WindowMasker](https://github.com/goeckslab/WindowMasker) to mark potential repeats on the genome. The genome is chunked into 10kb bins which move along the entire genome as sliding windows in order to profile the repeat intensity. Bedtools is then used to intersect the bins and WindowMasker fragments. These fragments are then mapped back to the original assembly for visualization purposes.

![Repeat Density workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_repeat_density.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### HIC_MAPPING

<details markdown="1">
<summary>Output files</summary>

- `hic_files/`
  - `*_pretext_hr.pretext`: High resolution pretext map.
  - `*_pretext_normal.pretext`: Low resolution pretext map.
  - `*.mcool`: HiC map required for HiGlass

</details>
The HIC_MAPPING subworkflow takes a set of HiC read files in .cram format as input and derives HiC mapping outputs in .pretext, .hic, and .mcool formats. These outputs are used for visualization on [PretextView](https://github.com/wtsi-hpag/PretextView), [Juicebox](https://github.com/aidenlab/Juicebox), and [Higlass](https://github.com/higlass/higlass) respectively.

![Hic Mapping workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_hic_mapping.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### TELO_FINDER

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`: A bgzipped file containing telomere sequence locations
  - `*.bed.gz.tbi`: A tabix index file for the above file.
- `hic_files/`
  - `*.bed`: The raw .bed file needed for ingestion into Pretext

</details>

The TELO_FINDER subworkflow uses a supplied (by the .yaml) telomeric sequence to identify putative telomeric regions in the input genome. This is acheived via the use of [FIND_TELOMERE_REGIONS](../modules/local/find_telomere_regions.nf), the output of which is used to generate a telomere.windows file with [FIND_TELOMERE_WINDOWS](../modules/local/find_telomere_windows.nf) (Both of these modules utilise VGP derived telomere programs [found here](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere)), data for each telomeric site is then extracted into bed format with [EXTRACT_TELO](../modules/local/extract_telo.nf) and finally BGZipped and indexed with [TABIX_BGZIPTABIX](https://nf-co.re/modules/tabix_bgziptabix/tabix_bgziptabix).

![Telomere Finder workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_telo_finder.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### BUSCO_ANALYSIS

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_buscogene.bigbed`: BigBed file for BUSCO genes track.
  - `*_ancestral.bigbed`: BigBed file for ancestral elements track.

</details>

The BUSCO_ANNOTATION subworkflow takes an assembly genome as input and extracts a list of [BUSCO](https://gitlab.com/ezlab/busco) genes based on the BUSCO results obtained from BUSCO. Additionally, it provides an overlap BUSCO gene set based on a list of lepidoptera ancestral genes (Wright et al., 2023), which has been investigated by Charlotte Wright from Mark Blaxter's lab at the Sanger Institute.

![Busco analysis workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_busco_analysis.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### GENE_ALIGNMENT

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

The gene alignment subworkflows loads genesets (cdna, cds, rna, pep) data from a given list of genomes detailed, in the input .yaml, and aligns these to the reference genome. It contains two subworkflows, one of which handles peptide data and the other of which handles RNA, nuclear and complementary DNA data. These produce files that can be displayed by JBrowse as tracks.

NUC_ALIGNMENTS: Reference fasta and fai files are aligned with the above mentioned gene alignment query files by [MINIMAP2_ALIGN](https://nf-co.re/modules/minimap2_align).
These are merged with [SAMTOOLS_MERGE](https://nf-co.re/modules/samtools_merge), converted to .bed format through [BEDTOOLS_BAMTOBED](https://nf-co.re/modules/bedtools_bamtobed), sorted via [BEDTOOLS_SORT](https://nf-co.re/modules/bedtools_sort) and finally converted to .bigBed format [UCSC_BEDTOBIGBED](https://nf-co.re/modules/ucsc_bedtobigbed) with the use of an auto SQL file found in the /assets/gene_alignment folder. This process is performed per species per data type.

PEP_ALIGNMENTS: Reference fasta is indexed with [MINIPROT_INDEX](https://nf-co.re/modules/miniprot_index) and aligned with peptide data [MINIPROT_ALIGN](https://nf-co.re/modules/miniprot_align). The output .gff file is merged with [CAT_CAT](https://nf-co.re/modules/cat_cat) per species, sorted with [BEDTOOLS_SORT](https://nf-co.re/modules/bedtools_sort) and indexed with [TABIX_BGZIPTABIX](https://nf-co.re/modules/tabix_bgziptabix/tabix_bgziptabix).

PUNCHLIST: Punchlists contain information on genes found to be duplicated (fully and partially) in the input genome. This is generated differently dependent on whether the datatype is peptide or not.

- NUC_ALIGNMENT:PUNCHLIST takes the merged.bam produced after the [SAMTOOLS_MERGE](https://nf-co.re/modules/samtools_merge) step. This is then converted into a .paf file with [PAFTOOLS_SAM2PAF](https://github.com/nf-core/modules/tree/master/modules/nf-core/paftools/sam2paf) and finally into bed with [PAF2BED](../modules/local/paf_to_bed.nf).
- PEP_ALIGNMENT:PUNCHLIST takes the merged.gff produced by [CAT_CAT](https://nf-co.re/modules/cat_cat) and converts it into .bed with [GFF_TO_BED](../modules/local/gff_to_bed.nf)

![Gene alignment workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_gene_alignment.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### INSILICO_DIGEST

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `{BSPQI|BSSSI|DLE1}.bigBed`: Bionano insilico digest cut sites track in the bigBed format for each of the set digestion enzymes.

</details>

The INSILICO_DIGEST workflow is used to visualize the Bionano enzyme cutting sites for a genome FASTA file. It starts by identifying the recognition sequences of the labeling enzyme to create a CMAP file. This CMAP file is then converted into BED and BIGBED formats to provide visualizations of the Bionano enzyme cutting sites. This procedure generates data tracks based on three digestion enzymes: BSPQ1, BSSS1, and DLE1.

![Insilico digest workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_insilico_digest.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### SELFCOMP

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*_selfcomp.bigBed`: BigBed file containing selfcomp track data.

</details>

he SELFCOMP subworkflow is a comparative genomics analysis originally performed by the Ensembl project. It involves comparing the genes and genomic sequences within a single species. The goal of the analysis is mainly to identify haplotypic duplications in a particular genome assembly.

![Selfcomp workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_selfcomp.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### SYNTENY

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.paf`: .paf file for each syntenic genomic aligned to reference.

</details>

This worflows searches along predetermined path for syntenic genome files based on clade and then aligns with [MINIMAP2_ALIGN](https://nf-co.re/modules/minimap2_align) each to the reference genome, emitting an aligned .paf file for each.

![Synteny workflow](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_synteny.jpeg)

![Workflow Legend](https://raw.githubusercontent.com/sanger-tol/treeval/dev/docs/images/treeval_1_0_legend.jpeg)

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
