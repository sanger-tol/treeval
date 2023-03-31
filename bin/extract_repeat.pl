#!/usr/local/bin/perl

use strict;

my $file = shift;
open(IN,$file);
my $last;
while (<IN>) {
    if (/\>(\S+)/) {
        $last = $1;
    }
    elsif (/(\d+)\s+-\s+(\d+)/) {
        print "$last\t$1\t$2\n";
    }
    else {
        die("Eh? $_");
    }
}
