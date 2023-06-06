# nf-core/treeval: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following workflows:

- [YAML_INPUT](#yamlinput) - Reads the input yaml and generates parameters used by other workflows.
- [GENERATE_GENOME](#generategenome) - Builds genome description file of the reference genome.
- [LONGREAD_COVERAGE](#longreadcoverage) - .
- [GAP_FINDER](#gapfinder) - .
- [REPEAT_DENSITY](#repeatdensity) - .
- [HIC_MAPPING](#hicmapping) - Aligns illumina HiC short reads to the input genome, generates mapping file in three format for visualisation: .pretext, .hic and .mcool
- [TELO_FINDER](#telofinder) - .
- [GENE_ALIGNMENT](#genealignment) - Aligns the peptide and nuclear data from assemblies of related species to the input genome.
- [INSILICO_DIGEST](#insilicodigest) - Generates a map of enzymatic digests using 3 Bionano enzymes.
- [SELFCOMP](#selfcomp) - Identifies regions of self-complementary sequence.
- [SYNTENY](#synteny) - Generates syntenic alignments between other high quality genomes.
- [BUSCO_ANALYSIS](#buscoanalysis) - Uses BUSCO to identify ancestral elements. Also use to identify ancestral Lepidopteran genes (merian units).

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### YAML_INPUT

This subworkflow reads the input .yaml via the use of the built-in snakeyaml.Yaml component, which converts the yaml into a nested list. Via some simple channel manipulation, each item in this nexted list is converted into a parameter for use in each of the other subworkflows.

### GENERATE_GENOME

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `my.genome`: Genome description file of the reference genome.

</details>

This workflow generates a .genome file which describes the base pair length of each scaffold in the reference genome. This is performed by [SAMTOOLS_FAIDX](https://nf-co.re/modules/samtools_faidx) to generate a .fai file. This index file is trimmed using local module [GENERATE_GENOME_FILE](../modules/local/generate_genome_file.nf) to output a .genome file. This file is then recycled into the workflow to be used by a number of other subworkflows.

<!--TODO: UPDATE FILE-->
![Generate genome workflow](images/treeval_generategenome_workflow.jpeg)


### LONGREAD_COVERAGE

<details markdown="1">
<summary>Output files</summary>

  - 

</details>


### GAP_FINDER

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`:
  - `*.bed.gz.tbi`:
- `hic_files/`
  - `*.bed`: The raw bed file needed for ingestion into Pretext

</details>

The GAP_FINDER subworkflow generates a bed file containing the genomic locations of the gaps in the sequence. This is performed by the use of [SEQTK_CUTN]() which cuts the input genome at sites of N (gaps). [GAP_LENGTH]() then calculates the lengths of gaps generates in the previous step, this file is injected into the hic_maps at a later stage. SEQTK's output bed file is then BGzipped and indexed by [TABIX_BGZIPTABIX](https://nf-co.re/modules/tabix_bgziptabix/tabix_bgziptabix).

<!-- ADD IMAGE -->

### REPEAT_DENSITY

<details markdown="1">
<summary>Output files</summary>

  - `hic_files/`
    - `coverage.bigWig`

</details>


### HIC_MAPPING

<details markdown="1">
<summary>Output files</summary>

  - `hic_files/`

</details>
The HIC_MAPPING subworkflow takes a set of HiC read files in CRAM format as input and derives HiC mapping outputs in .pretext, .hic, and .mcool formats. These outputs are used for visualization on PretextView (https://github.com/wtsi-hpag/PretextView), Juicebox (https://github.com/aidenlab/Juicebox), and Higlass (https://github.com/higlass/higlass) respectively.
The main steps involved include:

1. BWAMEM2_INDEX: This step indexes the input data using BWAMEM2. The output is redirected to a folder with the prefix BWAMEM2, which serves as a parameter for the mapping process.

2. CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT: This step is a complex process aimed at optimizing the performance of bwa-mem2 mem. It processes 10,000 containers from input CRAM files at a time and excludes the 5' chimeric reads. The mapping results also go through samtools fixmate to fill in information (insert size, cigar, mapq) about paired-end reads onto their corresponding other read. The final output is in BAM files.

3. The mapped BAM files are merged using SAMTOOLS_MERGE and fed into downstream processes:

4. PRETEXTMAP: This process generates pretext files based on the merged BAM files.

5. SAMTOOLS_MARKDUP: This process marks duplicate alignments in the merged BAM file.

6. The duplicate-marked BAM file is then converted to BED format and sorted using BAMTOBED_SORT.

7. Additionally, the paired contact reads are extracted using GET_PAIRED_CONTACT_BED based on the extracted paired contacts.

8. Finally, the extracted contacts are used to generate .hic and .mcool files using JUICER_TOOLS_PRE and [COOLER] respectively.
9. 
### TELO_FINDER

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`
  - `*.bed.gz`: A bgzipped file containing telomere sequence locations
  - `*.bed.gz.tbi`: A tabix index file for the above file.
- `hic_files/`
  - `*.bed`: The raw bed file needed for ingestion into Pretext

</details>

The TELO_FINDER subworkflow uses a suplied (by the .yaml) telomeric sequence to indentify putative telomeric regions in the input genome. This is acheived via the use of [FIND_TELOMERE_REGIONS](../modules/local/find_telomere_regions.nf), the output of which is used to generate a telomere.windows file with [FIND_TELOMERE_WINDOWS](../modules/local/find_telomere_windows.nf) (Both of these modules utilise VGP derived telomere programs [found here](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere)), data for each telomeric site is then extracted into bed format with [EXTRACT_TELO](../modules/local/extract_telo.nf) and finally BGZipped and indexed with [TABIX_BGZIPTABIX](https://nf-co.re/modules/tabix_bgziptabix/tabix_bgziptabix).

<!--ADD Figure-->


### BUSCO_ANALYSIS

<details markdown="1">
<summary>Output files</summary>

- `treeval_upload/`


</details>


### GENERATE_ALIGNMENT

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

<!--TODO: UPDATE FILE-->
![Gene alignment workflow](images/treeval_genealignment_workflow.jpeg)

### INSILICO_DIGEST

<details markdown="1">
<summary>Output files</summary>

- `insilico/`
  - `*.bigBed`: Bionano insilico digest cut sites track in the bigBed format for each of the set digestion enzymes.

</details>

This process runs for each of the digestion enzymes (bspq1, bsss1, DLE1). Using local module MAKECMAP_FA2CMAPMULTICOLOR to convert reference genome fasta into a colour-aware bionano .cmap format and emits files containing the index IDs and original genomic locations, which are passed into local module MAKECMAP_RENAMECMAPIDS to rename the .cmap IDs. This is used to create the .bed file (via MAKECMAP_CMAP2BED) and subsequently the .bigBed file (by [UCSC_BEDTOBIGBED](https://nf-co.re/modules/ucsc_bedtobigbed)) to be displayed as a JBrowse track.

<!--TODO: UPDATE FILE-->

![Insilico digest workflow](images/treeval_insilicodigest_workflow.jpeg)

### SELFCOMP

<details markdown="1">
<summary>Output files</summary>

- `selfcomp/`
  - `*.bigBed`: BigBed file containing selfcomp track data.

</details>

The reference fasta is split (SELFCOMP_SPLITFASTA) and chunked (CHUNKFASTA) to be by rapidly aligned with itself using [MUMMER](https://nf-co.re/modules/mummer). The outputted alignment files are merged (CONCATMUMMER) and converted into the .bed format (SELFCOMP_MUMMER2BED). This is then used by SELFCOMP_MAPIDS to generate a .bed file with a list of IDs and the genomic positions of selfcomplementary regions, which is then sorted by [BEDTOOLS_SORT](https://nf-co.re/modules/bedtools_sort). SELFCOMP_ALIGNMENTBLOCKS runs on this output to build alignment blocks. Merge alignment blocks ([BEDTOOLS_MERGE](https://nf-co.re/modules/bedtools_merge)) and then all individual block files (CONCATBLOCKS), filtered by motif length. This is converted to .bigBed by [UCSC_BEDTOBIGBED](https://nf-co.re/modules/ucsc_bedtobigbed).

<!--TODO: UPDATE FILE-->

![Selfcomp workflow](images/treeval_selfcomp_workflow.jpeg)

### SYNTENY

<details markdown="1">
<summary>Output files</summary>

- `synteny/`
  - `*.paf`: PAF file for each syntenic genomic aligned to reference.

</details>

This worflows searches along predetermined path for syntenic genome files based on clade and then aligns with [MINIMAP2_ALIGN](https://nf-co.re/modules/minimap2_align) each to the reference genome, emitting an aligned .paf file for each.

<!--TODO: UPDATE FILE-->

![Synteny workflow](images/treeval_synteny_workflow.jpeg)


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
