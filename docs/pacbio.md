## PacBio Data

Before running the pipeline data has to be in the `fasta.gz` format. Because of the software we use this data with it must also be long-read data as well as single stranded. This means you could use ONT too (except duplex reads).

The below commands should help you convert from mapped bam to fasta.gz, or from fastq to fasta.

If your data isn't already in these formats, then let us know and we'll see how we can help.

### BAM -> FASTQ

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

### FASTQ -> FASTA

This command creates a `fasta` folder (to store our fasta files), moves into the `fastq` folder and then converts `fastq` to `fasta` using seqtk seq.

```bash
mkdir fasta
cd fastq

for i in *fq; do
  echo $i
  j=${i%.fq}
  echo $j
  seqtk seq -a $i > ../fasta/${j}.fasta
done
```

### FASTA -> FASTA.GZ

This simply gzips the fasta files.

```bash
for i in .fasta; do
  echo $i
  gzip $i
done
```

### Or if you're a command line ninja

```bash
samtools bam2fq {prefix}.bam| seqtk seq -a - | gzip - > {prefix}.fasta.gz
```
