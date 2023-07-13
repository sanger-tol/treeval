[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/sanger-tol/treeval)

## Introduction

**sanger-tol/treeval** is a bioinformatics best-practice analysis pipeline for the generation of data supplemental to the curation of reference quality genomes. This pipeline has been written to generate flat files compatible with [JBrowse2](https://jbrowse.org/jb2/).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

The treeval pipeline has a sister pipeline currently named [curationpretext](https://github.com/sanger-tol/curationpretext) which acts to regenerate the pretext maps and accessory files during genomic curation in order to confirm interventions. This pipeline is sufficiently different to the treeval implementation that it is written as it's own pipeline.
## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

Currently, it is advised to run the pipeline with docker or singularity as a couple of major modules do not currently have a conda env associated with them.

Now, you can run the pipeline using:

```bash
nextflow run main.nf -profile singularity --input treeval.yaml -entry {FULL|RAPID} --outdir {OUTDIR}
```

An example treeval.yaml can be found [here](assets/local_testing/nxOscDF5033.yaml).

Further documentation about the pipeline can be found in the following files: [usage](https://nf-co.re/treeval/usage), [parameters](https://nf-co.re/treeval/parameters) and [output](https://nf-co.re/treeval/output).

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

sanger-tol/treeval has been written by Damon-Lee Pointon (@DLBPointon), Yumi Sims (@yumisims) and William Eagles (@weaglesBio).

We thank the following people for their extensive assistance in the development of this pipeline:

<ul>
  <li>@gq1 - For building the infrastructure around TreeVal</li>
  <li>@ksenia-krasheninnikova - For help with C code implementation and YAML parsing</li>
  <li>@mcshane - For guidance on algorithms </li>
  <li>@muffato - For code reviews and code support</li>
  <li>@priyanka-surana - For help with the majority of code reviews and code support</li>
</ul>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!--TODO: Citation-->

If you use sanger-tol/treeval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX).

### Tools

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
