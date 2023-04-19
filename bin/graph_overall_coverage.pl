#!/usr/bin/env perl

# Script originally developed by Yumi Sims (yy5@sanger.ac.uk)

use warnings;

# my $file = shift;

my ($file) = @ARGV;

if (!@ARGV || ($ARGV[0] eq '--version')) {
    print "1.0\n";
    exit 0;
}

open (FILE, $file) || die "can't open file $file\n";

my %depthcount;
while (my $line = <FILE>) {
    chomp $line;
    my ($id, $start, $end, $depth) = split ("\t", $line);
    my $length = $end - $start;

    if ($depthcount{$depth}){
        $depthcount{$depth} += $length;
    }
    else {
        $depthcount{$depth} = $length;
    }
}

foreach my $depth (sort {$a<=>$b} keys %depthcount){
    print join("\t", $depth, $depthcount{$depth}) ."\n";
}
