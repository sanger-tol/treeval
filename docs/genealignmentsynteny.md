## Work through

Seeing as this can be quite a complicated to set up here's a walk through.

### Step 1 - Set up the directories

Lets set up the system as if we want to run it on a bird genome.

```bash

mkdir -p gene_alignment_prep/scripts/

cp treeval/bin/treeval-dataprep/* gene_alignment_prep/scripts/

mkdir -p gene_alignment_prep/raw_fasta/

mkdir -p gene_alignment_data/bird/csv_data/

mkdir -p synteny/bird/
```

The naming of the bird folder here is important, keep this in mind.

So now we have this structure:

```
~/treeval-resources
    │
    ├─ synteny/
    │   └─ bird/
    │
    ├─ gene_alignment_data/
    │   └─ bird/
    │         └─ csv_data/
    │
    └─ gene_alignment_prep/
        ├─ scripts/
        └─ raw_fasta/
```

### Step 2 - Download some data

First, let's download out sytenic alignment data. I think the zebra finch would be good against the Chicken.

```bash
cd  synteny/bird/

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/957/565/GCA_003957565.4_bTaeGut1.4.pri/GCA_003957565.4_bTaeGut1.4.pri_genomic.fna.gz -o bTaeGut1_4.fasta.gz

gunzip bTaeGut1_4.fasta.gz
```

This leaves us with a file called `bTaeGut1_4.fasta` the genomic assembly of bTaeGut1_4 (the Darwin Tree of Life ID for this species) also known as `taeniopygia guttata`, the Australian Zebra Finch.

Now lets move into the raw_data folder and download some data, this may take some time.

```bash
cd  ../../gene_alignment_prep/raw_data/

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_cds_from_genomic.fna.gz -o GallusGallus-GRCg7b.cds.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz -o GallusGallus-GRCg7b.cdna.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_protein.faa.gz -o GallusGallus-GRCg7b.pep.fasta.gz

curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_rna.fna.gz -o GallusGallus-GRCg7b.rna.fasta.gz

```

Now that's all downloaded we need to prep it. At this point it is all still gzipped (the .gz on the end denotes that the file is compressed) in this format we can't use it. So lets use some bash magic.

This is a for loop written in bash, this will look through the current folder we are in for files ending with .fasta.gz and then gunzip them. This unzips the file, uncompresses it so it is usable, and then runs our python script, `GA_data_prep.py`.

```bash
for i in *.fasta.gz; do
gunzip $i;
python3 GA_data_prep.py ${i/.gz} ncbi 10;
done
```

This python script command here, in english, means. Take the file, uncompressed, that I downloaded from ncbi and cut it into pieces. A fasta file looks something like this, with headers (lines starting with `>`) and sequence (The usual ATGC's):

```markdown
>SCAFFOLD_1_DATA_ON_METHODS
ATGCGCATGCATGCATCGACTTCGAGCATCGTAG
>SCAFFOLD_2_DATA_ON_METHODS
ACCAGTGCTAGCTAGCTACGTGTGGGTTTGCCCCGTTT
```
The headers here will be trimmed, to only essential data that you need in order to fine the sequence in your database of choice.

Fasta file may be made up of anywhere between 10's to many thousands of these. So in the case of our cdna and pep files they need to be cut up to let TreeVal have a chance in reading them all in a small time frame.

cds and rna files will be cut up into 1000 header, sequence pairs per file.
pep and cdna will be cut up by a number you give or by 100.

This is because the size of pep and cdna files are so much larger.

The smaller the number you chunk a file, the smaller the files you produce which means you will also make many more files so there is a trade off.

So will produce a large amount of output in your terminal. Looking like:

```bash
python3 ../scripts/GA_data_prep.py GallusGallus-GRCg7b.cds.fasta ncbi 100

Your using at least Version 3.6, You are good to go...
os imported
argparse imported
regex imported
WORKING ON:             cds--GallusGallus-GRCg7b
Records per file:       1000
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

This is pretty much telling us that, yes you have given me a file and for every 1000 (i'm ignoring the number you gave me because this isn't a pep or cdna file) header, sequence pairs I have come across I have made a new file found here. You'll notice that it has also generated a new set of folders. This is based off of how we have named the file.

If you now type `ls` you should see the files we have downloaded (`GallusGallus-GRCg7b.cds.fasta`) and the folder `GallusGallus`. This folder we can now move to its permanent home.

```bash
mv GallusGallus/ ../../gene_alignment_data/bird/
```

### Step 3 -- Generate the csv

This file will act as an index of all files we have produced in the gene_alignment_data folder, and thankfully is a very simple step.

```bash
cd ../
python3 scripts/GA_csv_gen.py /path/to/gene_alignment_data/
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

So what is happening is that it is moving through the directory, identifying each unique folder and generating a csv summarising the data found in those directories into a csv with the following information, this looks like:

```bash
head -n 5 /gene_alignment_data/bird/csv_data/Gallus_gallus.GRCg6a-data.csv

org,type,data_file
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus9002cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus28453cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus18005cds.MOD.fa
Gallus_gallus.GRCg6a,cds,/gene_alignment_data/bird/Gallus_gallus/Gallus_gallus.GRCg6a/cds/Gallus_gallus6001cds.MOD.fa
```

This is all useful for the pipeline which generates job ids based on the org column, groups files by org and type columns and then pulls data from the data_file.

### Step 4 -- Understand where we are at

So we have now generated the directory structure for gene_alignment_data. So now lets use what we know to fill out the yaml.

The yaml is a file that we need in order to tell the pipeline where everything is, an example can be found here: [EXAMPLE YAML](../assets/local_testing/nxOscDF5033.yaml).

Here we can see a number of fields that need to be filled out, the easiest being `synteny_genome_path` and `data_dir`. These refer to the directories we made earlier so we can replace them as such:

```yaml

alignment:
  data_dir: /FULL/PATH/TO/treeval-resources/gene_alignment_data/

synteny_genome_path: /FULL/PATH/TO/treeval-resources/synteny

```

I said earlier that the the fact we called a folder `bird` was important, this is because it now becomes our `classT`:

```yaml
classT: bird
```

During the running of the pipeline, this is appended onto the end of `data_dir` and `synteny_genome_path` in order to find the correct files to use. So now all of the files inside `/FULL/PATH/TO/treeval-resources/synteny/bird/ ` will be used for syntenic alignments. Likewise with our genomic_alignment_data, TreeVal will turn this into `/FULL/PATH/TO/treeval-resources/gene_alignment_data/bird/` and then appends `csv_data`.

In Step 3, we generated some files which will be living in our `/FULL/PATH/TO/treeval-resources/gene_alignment_data/bird/csv_data/` folder and look like `GallusGallus.GRCg7b-data.csv`. These (minus the `-data.csv`) will be what we enter into the `geneset` field in the yaml. The common_name is a field we don't currently use.

```yaml
alignment:
  data_dir: /FULL/PATH/TO/treeval-resources/gene_alignment_data/
  common_name: "" # For future implementation (adding bee, wasp, ant etc)
  geneset: "GallusGallus.GRCg7b"
```

However, what is cool about this field is that you can add as many as you want. So say you have the genomic_alignment_data for the Finch saved as `TaeniopygiaGuttata.bTaeGut1_4`. The geneset field becomes: `geneset: "GallusGallus.GRCg7b,TaeniopygiaGuttata.bTaeGut1_4"`


https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
