# sanger-tol/treeval

[![GitHub Actions CI Status](https://github.com/sanger-tol/treeval/actions/workflows/ci.yml/badge.svg)](https://github.com/sanger-tol/treeval/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/sanger-tol/treeval/actions/workflows/linting.yml/badge.svg)](https://github.com/sanger-tol/treeval/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/sanger-tol/treeval)

## Introduction

**sanger-tol/treeval [1.2.0 - Ancient Destiny-]** is a bioinformatics best-practice analysis pipeline for the generation of data supplemental to the curation of reference quality genomes. This pipeline has been written to generate flat files compatible with [JBrowse2](https://jbrowse.org/jb2/) as well as HiC maps for use in Juicebox, PretextView and HiGlass.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

You can also set up and attempt to run the pipeline here: https://gitpod.io/#https://github.com/BGAcademy23/treeval-curation
This is a gitpod set up for BGA23 with a version of TreeVal, although for now gitpod will not run a nextflow pipeline die to issues with using singularity. We will be replacing this with an AWS instance soon.

The treeval pipeline has a sister pipeline currently named [curationpretext](https://github.com/sanger-tol/curationpretext) which acts to regenerate the pretext maps and accessory files during genomic curation in order to confirm interventions. This pipeline is sufficiently different to the treeval implementation that it is written as it's own pipeline.

1. Parse input yaml ( YAML_INPUT )
2. Generate my.genome file ( GENERATE_GENOME )
3. Generate insilico digests of the input assembly ( INSILICO_DIGEST )
4. Generate gene alignments with high quality data against the input assembly ( GENE_ALIGNMENT )
5. Generate a repeat density graph ( REPEAT_DENSITY )
6. Generate a gap track ( GAP_FINDER )
7. Generate a map of self complementary sequence ( SELFCOMP )
8. Generate syntenic alignments with a closely related high quality assembly ( SYNTENY )
9. Generate a coverage track using PacBio data ( LONGREAD_COVERAGE )
10. Generate HiC maps, pretext and higlass using HiC cram files ( HIC_MAPPING )
11. Generate a telomere track based on input motif ( TELO_FINDER )
12. Run Busco and convert results into bed format ( BUSCO_ANNOTATION )
13. Ancestral Busco linkage if available for clade ( BUSCO_ANNOTATION:ANCESTRAL_GENE )
14. Count KMERs with FastK and plot the spectra using MerquryFK ( KMER )
15. Generate a coverge track using KMER data ( KMER_READ_COVERAGE )

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Currently, it is advised to run the pipeline with docker or singularity as a small number of major modules do not currently have a conda env associated with them.

Now, you can run the pipeline using:

```bash
# For the FULL pipeline
nextflow run main.nf -profile singularity --input treeval.yaml --outdir {OUTDIR}

# For the RAPID subset
nextflow run main.nf -profile singularity --input treeval.yaml -entry RAPID --outdir {OUTDIR}
```

An example treeval.yaml can be found [here](assets/local_testing/nxOscDF5033.yaml).

Further documentation about the pipeline can be found in the following files: [usage](https://pipelines.tol.sanger.ac.uk/treeval/dev/usage), [parameters](https://pipelines.tol.sanger.ac.uk/treeval/dev/parameters) and [output](https://pipelines.tol.sanger.ac.uk/treeval/dev/output).

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

sanger-tol/treeval has been written by Damon-Lee Pointon (@DLBPointon), Yumi Sims (@yumisims) and William Eagles (@weaglesBio).

We thank the following people for their extensive assistance in the development of this pipeline:

<ul>
  <li>@gq1 - For building the infrastructure around TreeVal and helping with code review</li>
  <li>@ksenia-krasheninnikova - For help with C code implementation and YAML parsing</li>
  <li>@mcshane - For guidance on algorithms </li>
  <li>@muffato - For code reviews and code support</li>
  <li>@priyanka-surana - For help with the majority of code reviews and code support</li>
</ul>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use sanger-tol/treeval for your analysis, please cite it using the following doi: [10.5281/zenodo.10047653](https://doi.org/10.5281/zenodo.10047653).

### Tools

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
