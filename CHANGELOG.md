# sanger-tol/treeval: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

### Parameters

| Old Parameter | New Parameter |
| ------------- | ------------- |
| -             | --input       |

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own Biocontainer. This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Module                         | Old Version | New Versions              |
| ------------------------------ | ----------- | ------------------------- |
| bedtools                       | -           | 2.31.0                    |
| busco                          | -           | 5.4.3                     |
| bwa-mem2                       | -           | 2.2.1                     |
| cat                            | -           | 2.3.4                     |
| cooler                         | -           | 0.9.2                     |
| gnu-sort                       | -           | 8.25                      |
| minimap2 + samtools            | -           | 2.24 + 1.14               |
| miniprot                       | -           | 0.5                       |
| mummer                         | -           | 3.23                      |
| paftools (minimap2 + samtools) | -           | 2.24 + 1.14               |
| pretextmap + samtools          | -           | 0.1.9=h9f5acd7_1 + 1.16.1 |
| samtools                       | -           | 1.17                      |
| seqtk                          | -           | 1.3                       |
| tabix                          | -           | 1.11                      |
| ucsc                           | -           | 377                       |
| windowmasker (blast)           | -           | 2.13.0                    |

### Fixed

### Dependencies

### Deprecated
