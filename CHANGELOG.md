# sanger-tol/treeval: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.2] - Ancient Destiny (H2)- [2025-01-30]

Our 5th release for sanger-tol/treeval, correcting a software bug inside PretextGraph.

### Enhancements & Fixes

- Correction to the PRETEXT_GRAPH module, remade pretextgraph container with newest version 0.0.7.
- Change the way it takes read files, these should now be declared in the input yaml file. Details in the usage document.
- read_data can now include a fofn (file of file names) where each line contains one read file.
- read data is now checked for extension (this will look for fa(sta).gz or fofn---containing fa(asta).gz files ).
- Converted shell block modules into script block modules

### Software dependencies

| Module                | Old Version     | New Versions    |
| --------------------- | --------------- | --------------- |
| pretextmap + samtools | 0.0.2-c4 + 1.17 | 0.0.3-c1 + 1.17 |

## [1.2.1] - Ancient Destiny (H1)- [2025-01-22]

Our 4th release for sanger-tol/treeval, focusing on refining methods.

### Enhancements & Fixes

- Spelling mistake for the steps parameter: `--steps gene_alignment` rather than `--steps gene_alignments`
- Remove extra characters
- Correction to the PRETEXT_GRAPH module.

### Software dependencies

| Module                | Old Version  | New Versions    |
| --------------------- | ------------ | --------------- |
| pretextmap + samtools | 0.0.3 + 1.17 | 0.0.2-c4 + 1.17 |

## [1.2.0] - Ancient Destiny - [2024-11-15]

Our 3rd release for sanger-tol/treeval.

### Enhancements & Fixes

- Togglable subworkflows
- Adds a JBrowse Only workflow (this will lead to an update to the FULL workflow which can now call JBROWSE_ONLY and RAPID).
- Updates to containers (local modules) to remove Anaconda dependencies following policy changes.
- Updates to modules to remove Anaconda dependencies following policy changes
  - The majority of these updates only remove the `default` channel from the environment.yml
- CONDA warnings for modules which cannot use CONDA.
- Removable of a liberal use of spaces.
- reformat_intersect was previously not outputing version data.
- Adding arch specification to Pretext GitHub actions runner. Hopefully this will stop the spurious errors we see on there.
- Addition of steps into schema.
- Adds \*ktab as an output.
- Adds \*bin as an output for faster downsteam map generation.
- Updated singularity containers
- Added `--metaeuk` to BUSCO_BUSCO, default was causing pipeline errors on Actions -- Needs more investigation.
- Replaced Pyfasta split (depreciated 6 years ago) with Seqkit split which is frequently updated and very fast.
- Allocated resource review

### Parameters

| Old Parameter | New Parameter |
| ------------- | ------------- |
| -             | --steps       |

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own Biocontainer. This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Module                                 | Old Version      | New Versions      |
| -------------------------------------- | ---------------- | ----------------- |
| bamtobed_sort ( bedtools + samtools )  | 2.31.0 + 1.17    |                   |
| bedtools                               | 2.31.1           | -                 |
| busco                                  | 5.5.0            | -                 |
| bwa-mem2                               | 2.2.1            |                   |
| cat                                    | 2.3.4            |                   |
| chunk_fasta ( pyfasta )                | 0.5.2-1          | REMOVED           |
| cooler                                 | 0.9.2            |                   |
| cram_filter_align_bwamem2_fixmate_sort | -                |                   |
| ^ ( samtools + bwamem2 ) ^             | 1.17 + 2.2.1     |                   |
| coreutils                              | 9.1              |                   |
| fastk                                  | 1.0.1            |                   |
| gcc                                    | 10.4.0           |                   |
| find_telomere_windows ( java-jdk )     | 8.0.112          |                   |
| generate_cram_csv ( samtools )         | 1.17             |                   |
| gnu-sort                               | 8.25             | 9.3               |
| juicer_tools_pre ( java-jdk )          | 8.0.112          |                   |
| perl                                   | 5.26.2           |                   |
| merquryfk                              | 1.0.1            |                   |
| minimap2 + samtools                    | 2.24 + 1.14      |                   |
| minimap2_index                         | 2.24             | 2.28              |
| miniprot                               | 0.11--he4a0461_2 |                   |
| mummer                                 | 3.23             |                   |
| paftools ( minimap2 + samtools )       | 2.24 + 1.14      |                   |
| pretextmap + samtools                  | 0.0.2 + 1.17     | 0.0.3 + 1.17      |
| python                                 | 3.9              | -                 |
| - pandas                               | 1.5.2            | -                 |
| samtools                               | 1.18             | 1.21              |
| selfcomp_splitfasta ( perl-bioperl )   | 1.7.8-1          |                   |
| seqtk                                  | 1.4              |                   |
| seqkit                                 | ADDED            | 2.9.0--h9ee0642_0 |
| tabix                                  | 1.11             |                   |
| ucsc                                   | 377              | 447               |
| windowmasker (blast)                   | 2.14.0           | 2.15.0            |

- busco is currently pinned to v5.5.0 - Upgrading v5.7.1 would cause github actions to crash. Further investigation needed.

## [1.1.1] - Ancient Aurora (H1) - [2024-04-26]

### Enhancements & Fixes

- Generate CRAM CSV fix to allow for multi-readgroup cram files
- Removing KMER_READCOV
- tmp directory was being used
- Output file adjustment (names and location)

## [1.1.0] - Ancient Aurora - [2024-04-26]

The second release for sanger-tol, created with the [nf-core](https://nf-co.re/) template.

This builds on the initial release by adding subworkflows which generate kmer based coverage tracks and a kmer spectra graph. There are also a number of updates to logic used throughout the pipeline, as well as to the resources required by a significant number of modules.

### Enhancements & Fixes

- Updates to the resource allocation methods used by a number of modules in the base.config.
- Added a flag to stop the usage of Juicer.
- Subworkflow to generate a kmer based coverage track.
- Subworkflow to generate/update a kmer spectra graph.
- Subworkflow to use minimap2 for HiC mapping, if selected.
- Subworkflow to use BWAmem2 for HiC mapping, if selected.
- Subworkflow to ingest Pretext accessory files into the Pretext file, simplifying post-TreeVal data manipulation.
- Updated the logic in use throughout the pipeline.
- Updated the modules.config to include some of the logic, cleaning the code.
- Updated the HiC subworkflow to include subsampling the HiC data for Juicer due to resource requirements with large amounts of data.
- Updated the YAML_INPUT subworkflow, this now contains "flags" to change some software options.
- Updated the data names in the input YAML to reduce confusion.
- Updated software (Pretext{View,Snapshot,Graph}) to allow for use on large genomes with big data.
  - Added associated patch files and cpu architecture files.
- Updated the minimap2 align module to remove samtools view in preference of paftools for our usecase.
- Updated the test.yml inline with the above changes.
- Updated the SELFCOMP subworkflow to allow for the parallelisation of the work on large genomes.
- Updated the READ_COVERAGE subworkflow to produce the scaffold based AVG coverage and STND coverage
- Updated Modules from NF-Core - mostly relates to module structure rather than software.
- Updated the SummaryStats output to include HiC container counts.
- Added -T / -t flags where possible to minimise the use of the /tmp directory.
- Replaced CONCAT_MUMMER with CATCAT for simplicity.
- Removed JUICER from the RAPID entrypoint.
- Removed the csi or tbi logic. CSI is now used by default, this simplified the workflow and enlarges the capacity to handle much larger genomes. The logic block previously required was then moved.
- Added NF-DOWNLOAD to the CI-CD due to an error that causes incomplete downloaded when downloading a number of images at the same time.
- Added the RAPID_TOL entry point which is more geared towards the requirements of Sanger.
- Fix a bug in build_alignment_blocks.py to avoid indexing errors happening in large genomes.
- Change output BEDGRAPH from EXTRACT_TELO module.

#### Hot Fix 1

- Generate CRAM CSV fix to allow for multi-readgroup cram files
- Removing KMER_READCOV
- tmp directory was being used
- Output file adjustment (names and location)

### Parameters

| Old Parameter | New Parameter |
| ------------- | ------------- |
| -             | --juicer      |

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own Biocontainer. This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Module                                 | Old Version  | New Versions     |
| -------------------------------------- | ------------ | ---------------- |
| bamtobed_sort ( bedtools + samtools )  | -            | 2.31.0 + 1.17    |
| bedtools                               | 2.31.0       | 2.31.1           |
| busco                                  | 5.4.3        | 5.5.0            |
| bwa-mem2                               | -            | 2.2.1            |
| cat                                    | -            | 2.3.4            |
| chunk_fasta ( pyfasta )                | -            | 0.5.2-1          |
| cooler                                 | -            | 0.9.2            |
| cram_filter_align_bwamem2_fixmate_sort | -            |                  |
| ^ ( samtools + bwamem2 ) ^             | -            | 1.17 + 2.2.1     |
| coreutils                              | -            | 9.1              |
| fastk                                  | -            | 1.0.1            |
| gcc                                    | 7.1.0        | 10.4.0           |
| find_telomere_windows ( java-jdk )     | -            | 8.0.112          |
| generate_cram_csv ( samtools )         | -            | 1.17             |
| gnu-sort                               | -            | 8.25             |
| juicer_tools_pre ( java-jdk )          | -            | 8.0.112          |
| perl                                   | -            | 5.26.2           |
| merquryfk                              | -            | 1.0.1            |
| minimap2 + samtools                    | -            | 2.24 + 1.14      |
| miniprot                               | -            | 0.11--he4a0461_2 |
| mummer                                 | -            | 3.23             |
| paftools ( minimap2 + samtools )       | -            | 2.24 + 1.14      |
| pretextmap + samtools                  | 0.1.9 + 1.17 | 0.0.2 + 1.17     |
| python                                 | 3.9          | -                |
| - pandas                               | 1.5.2        | -                |
| samtools                               | 1.17         | 1.18             |
| selfcomp_splitfasta ( perl-bioperl )   | -            | 1.7.8-1          |
| seqtk                                  | -            | 1.4              |
| tabix                                  | -            | 1.11             |
| ucsc                                   | -            | 377              |
| windowmasker (blast)                   | -            | 2.14.0           |

### Fixed

- Resource allocations being calculated incorrectly.
- Pretext bugs related to large data.

### Dependencies

### Deprecated

## [1.0.0] - Ancient Atlantis - [2023-06-27]

Initial release of sanger-tol/treeval, created with the [nf-core](https://nf-co.re/) template.

The essential pathways of the gEVAL pipeline have now been converted to Nextflow DSL2 from vr-runner, snakemake and wr. Of the original pipeline there is only Bionano left to implement.

### Enhancements & Fixes

- Updated to nf-core/tools template v2.8.0.
- Subworkflow to generate channels from input yaml.
- Subworkflow to generate genome summary file using samtools
- Subworkflow to generate busco gene tracks and ancestral busco mapping.
- Subworkflow to generate HiC maps with cooler, juicebox and pretext.
- Subworkflow to generate gene alignments using miniprot and minimap2.
- Subworkflow to generate insilico digest tracks.
- Subworkflow to generate longread coverage tracks from pacbio data.
- Subworkflow to generate punchlists detailing regions of interest in the genome.
- Subworkflow to generate repeat density tracks.
- Subworkflow to generate tracks detailing self complementary regions.
- Subworkflow to generate syntenic alignments to high quality genomes.
- Subworkflow to generate tracks containing telomeric sites.
- Custom Groovy for reporting to provide file metrics and resource usage.
- Citations and all docs (including walkthroughs).
- Added gitpod.yml for running in the cloud. This is the tutorial written for BGA23.

### Parameters

| Old Parameter | New Parameter |
| ------------- | ------------- |
| -             | --input       |

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own Biocontainer. This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Module                                 | Old Version | New Versions     |
| -------------------------------------- | ----------- | ---------------- |
| assign_ancestal ( pandas + Python )    | -           | 1.5.2 + 3.9      |
| bamtobed_sort ( bedtools + samtools )  | -           | 2.31.0 + 1.17    |
| bedtools                               | -           | 2.31.0           |
| busco                                  | -           | 5.4.3            |
| bwa-mem2                               | -           | 2.2.1            |
| cat                                    | -           | 2.3.4            |
| chunk_fasta ( pyfasta )                | -           | 0.5.2-1          |
| cooler                                 | -           | 0.9.2            |
| concat_block ( coreutils )             | -           | 9.1              |
| concat_mummer ( coreutils )            | -           | 9.1              |
| cram_filter_align_bwamem2_fixmate_sort | -           |                  |
| ^ ( samtools + bwamem2 ) ^             | -           | 1.16.1 + 2.2.1   |
| extract_ancestral ( python )           | -           | 3.9              |
| extract_buscogene ( coreutils )        | -           | 9.1              |
| extract_cov_id ( coreutils )           | -           | 9.1              |
| extract_repeat ( perl )                | -           | 5.26.2           |
| extract_telo ( coreutils )             | -           | 9.1              |
| find_telomere_regions ( gcc )          | -           | 7.1.0            |
| find_telomere_windows ( java-jdk )     | -           | 8.0.112          |
| findhalfcoverage ( python )            | -           | 3.9              |
| gap_length ( coreutils )               | -           | 9.1              |
| generate_cram_csv ( samtools )         | -           | 1.17             |
| get_largest_scaff ( coreutils )        | -           | 9.1              |
| get_paired_contact_bed ( coreutils )   | -           | 9.1              |
| get_synteny_genomes ( coreutils )      | -           | 9.1              |
| getminmaxpunches ( coreutils )         | -           | 9.1              |
| graphoverallcoverage ( perl )          | -           | 5.26.2           |
| gnu-sort                               | -           | 8.25             |
| juicer_tools_pre ( java-jdk )          | -           | 8.0.112          |
| makecmap_cmap2bed ( python )           | -           | 3.9              |
| makecmap_fa2cmapmulticolor ( perl )    | -           | 5.26.2           |
| makecmap_renamecmapids ( perl )        | -           | 5.26.2           |
| minimap2 + samtools                    | -           | 2.24 + 1.14      |
| miniprot                               | -           | 0.11--he4a0461_2 |
| mummer                                 | -           | 3.23             |
| paf_to_bed ( coreutils )               | -           | 9.1              |
| paftools ( minimap2 + samtools )       | -           | 2.24 + 1.14      |
| pretextmap + samtools                  | -           | 0.1.9 + 1.17     |
| reformat_intersect ( coreutils )       | -           | 9.1              |
| reformat_ids ( coreutils )             | -           | 9.1              |
| replace_dots ( coreutils )             | -           | 9.1              |
| samtools                               | -           | 1.17             |
| selfcomp_alignmentblocks ( python )    | -           | 3.9              |
| selfcomp_mapids ( python )             | -           | 3.9              |
| selfcomp_mummer2bed ( python )         | -           | 3.9              |
| selfcomp_splitfasta ( perl-bioperl )   | -           | 1.7.8-1          |
| seqtk                                  | -           | 1.4              |
| tabix                                  | -           | 1.11             |
| ucsc                                   | -           | 377              |
| windowmasker (blast)                   | -           | 2.14.0           |

### Fixed

### Dependencies

### Deprecated
