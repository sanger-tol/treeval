#!/usr/bin/env perl

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

#
# Take a genome fasta file, and creates a fasta file and agp
# that is suitable for loading into Ensembl
#
# If the genome is already suitable for loading (i.e. all of the
# sequences are small enough for direct loading), then the script
# says so and does nothing
#
# author: grit@sanger.ac.uk
use strict;

use Bio::SeqIO;
use Getopt::Long;

my $MAX_CONTIG_LEN = 500000;

my ($genome_fa,
    $contig_fa_file,
    $agp_file) = @ARGV;

if (!@ARGV || ($ARGV[0] eq '--version')) {
    print "1.0\n";
    exit 0;
}

my $seqio = Bio::SeqIO->new(-format  => 'fasta',
                            -file => ($genome_fa =~ /\.gz$/) ? "gunzip -c $genome_fa |" : $genome_fa,
);

my (@toplevels, $need_agp, $count);

while (my $seq = $seqio->next_seq) {
    if ($seq->length > $MAX_CONTIG_LEN) {
        $need_agp = 1;
    }
    push @toplevels, $seq;
}

if (not $need_agp) {
    print "All sequences are small enough for direct loading. Did not create AGP\n";
} else {
    my $outseqio = Bio::SeqIO->new(-format => 'fasta',
                                    -file => ">$contig_fa_file");
    open(my $agp_fh, ">$agp_file") or die "Could not open $agp_file for writing\n";

    foreach my $seq (@toplevels) {
        if ($seq->length < $MAX_CONTIG_LEN) {
            $outseqio->write_seq($seq);
            printf($agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n", $seq->id, 1, $seq->length, ++$count, $seq->id, 1, $seq->length);
        } else {
            print STDERR "GOT ONE\n";
            my $seg_count = 1;
            for(my $i = 0; $i < $seq->length; $i += $MAX_CONTIG_LEN) {
                my $start = $i + 1;
                my $end = $start + $MAX_CONTIG_LEN - 1;
                $end = $seq->length if $end > $seq->length;

                my $new_id = sprintf("%s.%d", $seq->id, $seg_count++);

                my $seq_seg = Bio::PrimarySeq->new( -id => $new_id,
                                                    -seq => substr($seq->seq, $start - 1, $end - $start + 1));
                $outseqio->write_seq($seq_seg);

                printf($agp_fh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t+\n", $seq->id, $start, $end, ++$count, $seq_seg->id, 1, $seq_seg->length);
            }
        }
    }
}

exit(0);
