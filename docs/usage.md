# sanger-tol/treeval: Usage

## :warning: Please read this documentation on the sanger-tol website: [https://pipelines.tol.sanger.ac.uk/treeval/](https://pipelines.tol.sanger.ac.uk/treeval/)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The TreeVal pipeline has a few requirements before being able to run:

- The `gene_alignment_data` requires a specific .csv format.

- HiC CRAM files must be pre-indexed in the same location as the CRAM file, e.g., `samtools index {cram file}`. A check and automated indexing of the cram file will be added in the future.

- Finally, the yaml file (which is described below in Full Samplesheet). This needs to contain all of the information related to the assembly for the pipeline to run.

## Prior to running TreeVal

:warning: Please ensure you read the following sections on Directory Structure (`gene_alignment_data` and scripts), HiC data prep and Pacbio data prep. Without these you may not be able to successfully run the TreeVal pipeline. If nothing is clear then please leave an issue report.

We now also support ( and encourage ) using the nf-co2footprint plugin (on Nextflow versions >= 23.07) which generates statistics on how much energy your pipeline uses as well as the amount of Co2 it helps produce. As it is pre-release, you will need to compile this plugin your self and store it in your `$NXF_HOME/plugins` directory, which you can find with `echo $NXF_HOME`. We have included the relevant config file `co2footprint.config` in this repo. The plugin can be used be including `-plugins nf-co2footprint@{VERSION} -c co2footprint.config` in your nextflow command. Please head to the website to find out more [NF-CO2FOOTPRINT](https://nextflow-io.github.io/nf-co2footprint/contributing/setup/).

### Local testing

<details markdown="1">
  <summary>Details</summary>

We provide a complete set of test data that can be used to test the pipeline locally.

```bash
git clone https://github.com/sanger-tol/treeval.git
cd treeval
curl https://tolit.cog.sanger.ac.uk/test-data/resources/treeval/TreeValTinyData.tar.gz | tar xzf -

sed -i "s|/home/runner/work/treeval/treeval|${PWD}|" TreeValTinyData/gene_alignment_data/fungi/csv_data/LaetiporusSulphureus.gfLaeSulp1-data.csv
sed -i "s|/home/runner/work/treeval/treeval|${PWD}|" assets/github_testing/TreeValTinyFullTest.yaml
```

This downloads the repo and test-data, which is then de-compressed. The sed commands then rewrites references from GitHub locations to local locations in the gene_alignment csv file as well as the treeval yaml file. The above command

You should now be able to run the pipeline as you see fit, something like the below makes sense in this case

```bash
nextflow run main.nf -profile test_github,singularity
```

</details>

### Gene Alignment and Synteny Data

<details markdown="1">
  <summary>Details</summary>

#### Step 1 -- Preparing Synteny data

For synteny you should provide the full genomic fasta file, of any high quality genome you want to be compared against.

For bird we recommend the Golden Eagle ( _Aquila chrysaetos_ ) and the Zebrafinch (_Taeniopygia guttata_), which can be downloaded from NCBI.

First, lets quickly make some folders

```bash
mkdir -p synteny/bird/
mkdir -p gene_alignment_prep/raw_data/
mkdir -p gene_alignment_data/bird/
```

Now, let's download some syntenic alignment data. I think the Zebrafinch (_Taeniopygia guttata_) would be good example.

```bash
cd  synteny/bird/

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.4_bTaeGut1.4.pri/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz -o bTaeGut1_4.fasta.gz

gunzip bTaeGut1_4.fasta.gz
```

This leaves us with a file called `bTaeGut1_4.fasta` the genomic assembly of `bTaeGut1_4` (this is a ToLID which you can read more about here: [Tree of Life ID](https://id.tol.sanger.ac.uk)) also known as _Taeniopygia guttata_, the Australian Zebrafinch.

Now lets move into the `/gene_alignment_prep/raw_data/` folder and download some data, this may take some time.

#### Step 1 -- Preparing Gene alignment data

```
cd  ../../gene_alignment_prep/raw_data/

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_cds_from_genomic.fna.gz -o GallusGallus-GRCg7b.cds.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz -o GallusGallus-GRCg7b.cdna.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa.gz -o GallusGallus-GRCg7b.pep.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_rna.fna.gz -o GallusGallus-GRCg7b.rna.fasta.gz
```

Keep note of how these files have been downloaded, the scripts used here need the file name to be in this format: `GallusGallus-GRCg7b.rna.fasta.gz` or `{Organism}-{Accession}.{data-type}.fasta.gz`

Now that's all downloaded we need to prep it. At this point it is all still gzipped (the `.gz` on the end denotes that the file is compressed) in this format we can't use it.

The below code will look through the current folder for files ending with `.fasta.gz` and decompresses it, it will then run our python script, `GA_data_prep.py`.

NOTE: This will soon be superseeded by a rust tool (`treeval_utils`), this will replace all of the accessory pythons scripts included in this document. Whilst also adding a program to verify your YAML.

Command:

```bash
for i in *.fasta.gz; do
  gunzip $i;
  python3 GA_data_prep.py ${i/.gz} ncbi 10;
done
```

The Python script will clean the headers common in NCBI and ensembl fasta files and then split the files into 10, or however many you choose, header/sequence pairs per file.

A fasta file may be made up of anywhere between 10's to many thousands of these pairs. So in the case of our `cdna` and `pep` files they need to be cut up to let TreeVal have a chance in reading them all in a small time frame. `cds` and `rna` files will be cut up into 1,000 header-sequence pairs per file. The number given on the command-line is ignored. `pep` and `cdna` will be cut up by a number you give or by 100. This is because the size of `pep` and `cdna` files are so much larger.

The smaller the number you chunk a file, the smaller the files you produce which means you will also make many more files so there is a trade off.

The `GA_data_prep.py` will produce a large amount of output in your terminal. Looking like:

```bash
python3 ../../bin/treeval-dataprep/GA_data_prep.py GallusGallus-GRCg7b.cds.fasta ncbi 100

Your using at least Version 3.6, You are good to go...
os imported
argparse imported
regex imported
WORKING ON: cds--GallusGallus-GRCg7b
Records per file: 100
Entryfunction called
GallusGallus-GRCg7b.cds.fasta
File found at %s GallusGallus-GRCg7b.cds.fasta
Renaming headers
Read_fasta called
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus1000cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus2001cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus3002cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus4003cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus5004cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus6005cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus7006cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus8007cds.MOD.fa
File saved: -- ./GallusGallus/GallusGallus.GRCg7b/cds/GallusGallus9008cds.MOD.fa
```

This is pretty much telling us that, yes you have given me a file and for every 100 (i'm ignoring the number you gave me because this isn't a `pep` or `cdna` file) header, sequence pairs I have come across I have made a new file found here. You'll notice that it has also generated a new set of folders. This is based off of how we have named the file.

If you now type `ls` you should see the files we have downloaded (`GallusGallus-GRCg7b.cds.fasta`) and the folder `GallusGallus`. This folder can now be moved to its permanent home.

```bash
mv GallusGallus/ ../bird/
```

Yes, we could have made this folder right away inside the `/bird` folder, however, it is nice to have that split of unprocessed and processed reads in case something unexpected happens.

#### Step 3 -- Generate the CSV

This file will act as an index of all files we have produced in the `gene_alignment_data` folder, and thankfully is a very simple step.

```bash
cd ../../ # So that we are now in the main treeval folder
python3 bin/treeval-dataprep/GA_csv_gen.py /gene_alignment_data/
```

Running this will look like:

```bash
============> CorvusMon1.bCorMon1 -- bird
Generating CSV for:     CorvusMon1.bCorMon1
Save Path:              /gene_alignment_data/bird/csv_data/CorvusMon1.bCorMon1-data.csv
============> CorvusMoneduloides.bCorMon1 -- bird
Generating CSV for:     CorvusMoneduloides.bCorMon1
Save Path:              /gene_alignment_data/bird/csv_data/CorvusMoneduloides.bCorMon1-data.csv
============> Gallus_gallus.UW_022020 -- bird
Generating CSV for:     Gallus_gallus.UW_022020
Save Path:              /gene_alignment_data/bird/csv_data/Gallus_gallus.UW_022020-data.csv
============> Gallus_gallus.GRCg6a -- bird
Generating CSV for:     Gallus_gallus.GRCg6a
Save Path:              /gene_alignment_data/bird/csv_data/Gallus_gallus.GRCg6a-data.csv
============> GallusGallus.GRCg7b -- bird
Generating CSV for:     GallusGallus.GRCg7b
Save Path:              /gene_alignment_data/bird/csv_data/GallusGallus.GRCg7b-data.csv
```

So what is happening is that it is moving through the directory tree of the gene_alignment folder, identifying each unique folder and generating a CSV summarising the data found in those directories into a csv with the following information:

```bash
head -n 5 /gene_alignment_data/bird/csv_data/Gallus_gallus.GRCg6a-data.csv

org,type,data_file
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus9002cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus28453cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus18005cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus6001cds.MOD.fa
```

This is all useful for the pipeline which generates job ids based on the org column, groups files by org and type columns and then pulls data from the data file.

#### Step 4 -- Understand where we are at

Now let's use what we know to fill out the yaml.

The yaml is a file that we need in order to tell the pipeline where everything is, an example can be found [here](https://raw.githubusercontent.com/sanger-tol/treeval/dev/assets/local_testing/nxOscDF5033.yaml).

```yaml
alignment:
  genesets:
    - /FULL/PATH/TO/<geneset_name>-data.csv
synteny:
  - /FULL/PATH/TO/<genome_name>.fasta
```

</details>

### HiC data Preparation

<details markdown="1">
  <summary>Details</summary>

Illumina HiC read files should be presented in an unmapped CRAM format, each must be accompanied by an index file (.crai) generated by samtools index. If your unmapped HiC reads are in FASTQ format, you should first convert them to CRAM format by using samtools import methods. Examples are below:

#### Conversion of FASTQ to CRAM

```bash
samtools import -@8 -r ID:{prefix} -r CN:{hic-kit} -r PU:{prefix} -r SM:{sample_name} {prefix}_R1.fastq.gz {prefix}_R2.fastq.gz -o {prefix}.cram
```

#### Indexing of CRAM

```bash
samtools index {prefix}.cram
```

</details>

### Longread Data Preparation

<details markdown="1">
  <summary>Details</summary>

Before running the pipeline, longread data must to be in the `fasta.gz` format. Because of the software we use this data with, it must also be long-read data and single stranded. This means you can use ONT too (except duplex reads), tested as of Feb/2025.

The below commands should help you convert from mapped bam to fasta.gz, or from fastq to fasta.

If your data isn't already in these formats, then let us know and we'll see how we can help.

#### BAM -> FASTQ

This command iterates through your bam files and converts them to fastq via samtools.

```bash
cd { TO FOLDER OF BAM FILES }
mkdir fastq
for i in *bam
do
  echo $i
  j=${i%.bam}
  echo $j
  samtools bam2fq ${i} > fastq/${j}.fq
done
```

#### FASTQ -> FASTA.GZ

This command creates a `fasta` folder (to store our fasta files), moves into the `fastq` folder and then converts `fastq` to `fasta` using `seqtk seq`.

```bash
mkdir fasta
cd fastq
for i in *fq; do
  echo $i
  j=${i%.fq}
  echo $j
  seqtk seq -a $i > ../fasta/${j}.fasta
  gzip ../fasta/${j}.fasta
done
```

#### Or if you're a command line ninja

You can do it all in one line, bam -> fasta.gz

```bash
samtools bam2fq {prefix}.bam | seqtk seq -a - | gzip - > {prefix}.fasta.gz
```

</details>

### Pretext Accessory File Ingestion

<details markdown="1">
  <summary>Details</summary>

Note: This will require you to install bigwigToBedGraph from the ucsc package. Instructions on downloading this can be found at [EXAMPLE #3](https://genome.ucsc.edu/goldenPath/help/bigWig.html#:~:text=Alternatively%2C%20bigWig%20files%20can%20be,to%20the%20Genome%20Browser%20server.)

The PreText files generated by the pipeline _are_ automatically ingested into the pretext files. However, you may want to ingest from other versions of pretextgraph or have your own pre-generated file you want to ingest. For this you must use the following code:

```
cd {outdir}/hic_files

bigWigToBedGraph {coverage.bigWig} /dev/stdout | PretextGraph -i { your.pretext } -n "coverage"

bigWigToBedGraph {repeat_density.bigWig} /dev/stdout | PretextGraph -i { your.pretext } -n "repeat_density"

cat {telomere.bedgraph} | awk -v OFS="\t" '{$4 = 1000; print}'|PretextGraph -i { your.pretext } -n "telomere"

cat {gap.bedgraph} | awk -v OFS="\t" '{$4= 1000; print}'| PretextGraph -i { your.pretext } -n "gap"
```

</details>

## Full samplesheet

YAML is "Yet Another Markdown Language", it is a human-readable format that we use to tell TreeVal a number of things. This includes; assembly location, telomere motif, longread data files (in fasta.gz format) and HiC cram files. The full Yaml is detailed below.

### YAML contents

The following is an example YAML file we have used during production: [nxOscDF5033.yaml](../assets/local_testing/nxOscDF5033.yaml) and is shown below. This contains some annotations we believe to be helpful, information on the alignment, synteny, longread and hic data.

- `assembly`
  - `assem_level`: scaffold or contig level assembly (not used).
  - `assem_version`: Used to complete sample_id.
  - `sample_id`: ToLID of the sample.
  - `latin_name`: Latin identification of species
  - `defined_class`: Clade name (as used to group synteny sequences and to complete alignment/data_dir).
  - `project_id`: Project id for the ticket (not used)
- `reference_file`: Sample .fa file.
- `assem_reads`
  - `read_type`: { hifi | clr | ont | illumina } To be used in future update.
  - `read_data`:
    - List of paths (ending with `/`) to folder containing fasta.gz files.
  - `supplementary_data`: Will be required in future development.
- `hic_data`:
  - `hic_cram`: path (ending with `/`) to folder containing cram files.
  - `hic_aligner`: choice between `bwamam2` and `minimap2`
- `alignment`
  - `genesets`:
    - List of Gene alignment data .csv file paths.
- `kmer_profile`:
  - `kmer_length`: length of kmer to be used in plotting, normally 31
  - `dir`: directory containing old plot to be regenerated if applicable
- `self_comp`
  - `motif_len`: Length of motif to be used in self complementary sequence finding
- `synteny`
  - List of paths to syntenic genomes grouped by clade.
- `intron:`
  - `size`: base pair size of introns default is 50k
- `telomere`:
  - `teloseq`: Telomeric motif
- `busco`
  - `lineages_path`: path to folder above lineages folder
  - `lineage`: Example is `nematode_odb10`

<details markdown="1">
  <summary>Notes on using BUSCO</summary>

The pipeline requires the use of BUSCO odb database.
Create the database directory and move into the directory:

```bash
DATE=2025_02
BUSCO=/path/to/databases/busco_${DATE}
mkdir -p $BUSCO
cd $BUSCO
```

Download BUSCO data and lineages to allow BUSCO to run in offline mode.

```bash
wget -r -nH https://busco-data.ezlab.org/v5/data/
```

The trailing slash after data is important, otherwise wget doesn't download the subdirectories.

Tar gunzip (decompress) all folders that have been stored as tar.gz, in the same parent directories as where they were stored:

```bash
find v5/data -name "*.tar.gz" | while read -r TAR; do tar -C `dirname $TAR` -xzf $TAR; done
```

If you have [GNU parallel](https://www.gnu.org/software/parallel/) installed, you can also use the command below which will run faster as it will run the decompression commands in parallel:

```bash
find v5/data -name "*.tar.gz" | parallel "cd {//}; tar -xzf {/}"
```

</details>

## Sub-workflows

<details markdown="1">
  <summary>Sub-workflows</summary>

- `YAML_INPUT`
  - Reads the input yaml and generates parameters used by other workflows.
- `GENERATE_GENOME`
  - Builds genome description file of the reference genome.
- `READ_COVERAGE`
  - Produces read coverage based on pacbio long read fasta file.
- `GAP_FINDER`
  - Identifies contig gaps in the input genome.
- `REPEAT_DENSITY`
  - Reports the intensity of regional repeats within an input assembly.
- `HIC_MAPPING`
  - Aligns illumina HiC short reads to the input genome, generates mapping file in three format for visualisation: .pretext, .hic and .mcool
- `TELO_FINDER`
  - Find a user given motif in the input genome.
- `GENE_ALIGNMENT`
  - Aligns the peptide and nuclear data from assemblies of related species to the input genome.
- `INSILICO_DIGEST`
  - Generates a map of enzymatic digests using 3 Bionano enzymes.
- `SELFCOMP`
  - Identifies regions of self-complementary sequence.
- `SYNTENY`
  - Generates syntenic alignments between other high quality genomes.
- `BUSCO_ANALYSIS`
  - Uses BUSCO to identify ancestral elements. Also use to identify ancestral Lepidopteran genes (merian units).
- `KMER`
  - Generating kmer graphs of the assembly.
- `KMER_READ_COVERAGE`
  - Generating read coverage using kmers.

</details>

## Running the pipeline

The typical command for running the pipeline is as follows (if you require the RAPID workflow you can append `-entry RAPID` to the command):

```console
nextflow run sanger-tol/treeval --input assets/treeval.yaml --outdir <OUTDIR> -profile singularity,sanger
```

With the `treeval.yaml` containing the information from the above YAML Contents section.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull sanger-tol/treeval
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [sanger-tol/treeval releases page](https://github.com/sanger-tol/treeval/releases) and find the latest version number - numeric only (eg. `1.2.2`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.2.2`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/sanger-tol/blobtoolkit/blob/56906ffb5737e4b985797bb5fb4b9c94cfe69600/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this, test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
