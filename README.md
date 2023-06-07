[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/treeval)

## Introduction

**nf-core/treeval** is a bioinformatics best-practice analysis pipeline for the generation of data supplemental to the curation of reference quality genomes. This pipeline has been written to generate flat files compatable with [JBrowse2](https://jbrowse.org/jb2/).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

The version 1 pipeline will be made up of the following steps, (r) = Steps run in Rapid:

- INPUT_READ (r)
  > The reading of the input yaml and conversion into channels for the sub-workflows.

- GENERATE_GENOME (r)
  > Generate .genome for the input genome using SAMTOOLS FAIDX.

- GENERATE_ALIGNMENT
  > Peptides will run pep_alignment.nf with Miniprot.

  > CDNA, RNA and CDS will run through nuc_alignment.nf with Minimap2.

- INSILICO DIGEST
  > Generates a map of enzymatic digests using 3 Bionano enzymes.

- SELFCOMP
  > Identifies regions of self-complementary sequencs using Mummer.

- SYNTENY
  > Generates syntenic alignments between other high quality genomes via Minimap2.

- ANCESTRAL ELEMENT ANALYSIS
  > Lepidopteran Element Analysis. Using BUSCO and custom python scripts to parse ancestral lep genes. This will eventually have a number of clade specific sub-workflows.

- LONGREAD_COVERAGE (r)
  > Calculating the coverage of reads across the genome.

- FIND_GAPS (r)
  > Identifying gaps in the input genome using seqtk cutn.

- FIND_TELOMERE (r)
  > Identify sites of a given telomeric sequence.

- REPEAT_DENSITY (r)
  > Generate a graph showing the relative amount of repeat in a given chunk.

- HIC_MAPPING (r)
  > Generation of HiC maps for the curation of a genome, these include: pretext_hires, pretext_lowres and cooler maps.

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

Now, you can run the pipeline using:

```bash
nextflow run main.nf -profile singularity --input treeval.yaml -entry {FULL|RAPID} --outdir {OUTDIR}
```

## Documentation

The nf-core/treeval pipeline comes with documentation about the pipeline [usage](https://nf-co.re/treeval/usage), [parameters](https://nf-co.re/treeval/parameters) and [output](https://nf-co.re/treeval/output).

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

nf-core/treeval has been written by Damon-Lee Pointon (@DLBPointon), Yumi Sims (@yumisims) and William Eagles (@weaglesBio).

We thank the following people for their extensive assistance in the development of this pipeline:

<ul>
  <li>@muffato</li>
  <li>@gq1</li>
  <li>@ksenia-krasheninnikova</li>
  <li>@priyanka-surana</li>
  <li>@mcshane</li>
</ul>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!--TODO: Citation-->
If you use nf-core/treeval for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX).

### Tools

<ul>
  <li>BedTools</li>
  <li>BUSCO </li>
  <li>bwa-mem2</li>
  <li>Cooler</li>
  <li>Find Telomere<//li>
  <li>Juicer</li>
  <li>Minimap2</li>
  <li>Miniprot</li>
  <li>Mummer</li>
  <li>Python3</li>
  <ul>
    <li>Pyfasta</li>
    <li>Pandas</li>
    <li>Pybedtools</li>
  </ul>
  <li>Perl</li>
  <li>Pretext</li>
  <li>Samtools</li>
  <li>Staden io_lib</li>
  <li>Tabix</li>
  <li>UCSC_Tools - Kent</li>
  <li>WindowMasker</li>

</ul>

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
