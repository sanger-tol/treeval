#! /usr/local/bin/perl
use warnings;
use strict;

use Getopt::Long;

my ($cmapfile, $xmapfile, $smapfile, $indelfile, $namefile, $bedfile, $idx_key, $agp, $multi);

GetOptions ('cmapfile:s'    =>  \ $cmapfile,
            'xmapfile:s'    =>  \ $xmapfile,
            'smapfile:s'    =>  \ $smapfile,
            'indelfile:s'   =>  \ $indelfile,
            'namefile:s'    =>  \ $namefile,
            'bedfile:s'     =>  \ $bedfile,
            'idx_key:s'     =>  \ $idx_key,
            'agprev:s'      =>  \ $agp,
            'multi!'        =>  \ $multi, # rename xmap genome map ids for multiple mappings
    );

my $namedata;
if ($idx_key){
    if ($agp){
        $namedata = readfiles($idx_key, "agp");
    }
    else {
        $namedata =  readfiles($idx_key, "idx");
    }
}
else {
    $namedata = ($cmapfile) ? readfiles($namefile, "cmap") : readfiles($namefile, "xmap");
}

if ($cmapfile) {
    open (CMAP, $cmapfile);

    my %idx;
    foreach my $line (<CMAP>){
        chomp $line;

        if ($line =~ /#/){
            print $line."\n";
            next;
        }

        my @F = split (/\t/, $line);

        my $key;
        if ($idx_key){
            $key = $F[0];
        }
        else {
            ($key = $F[1]) =~ s/\.0//;
        }

        die "misssing name at $key\n" if !($$namedata{$key});
        #print $$namedata{$key}."-$key\n";

        print join ("\t", $$namedata{$key}, $F[1], $F[2], $F[3], $F[4], $F[5], $F[6], $F[7], $F[8])."\n";

        # -if no original key file avaialble)
        $idx{$F[0]} = { name   => $$namedata{$key},
                        length => $key,
        } if (!$idx{$F[0]} && !$idx_key);
    }
    close CMAP;

    # -if no original key file avaialble)
    if (!$idx_key){
        open (KEY, ">>key.out");

        foreach my $index (sort {$a <=> $b}  keys %idx){
            print KEY join ("\t", $index, $idx{$index}{name}, $idx{$index}{length}) ."\n";
        }
        close KEY;
    }
}
elsif ($xmapfile){

    open (XMAP, $xmapfile);

    my %idx;
    foreach my $line (<XMAP>){
        chomp $line;

        if ($line =~ /#/){
            print $line."\n";
            next;
        }
        my @F   =  split (/\t/, $line);

        print join ("\t", $F[0], $F[1], $$namedata{$F[2]}, $F[3], $F[4], $F[5], $F[6], $F[7], $F[8], $F[9], $F[10], $F[11], $F[12], $F[13])."\n";

    }
    close (XMAP);
}
elsif ($indelfile){

    open (INDEL, $indelfile);
    foreach my $line (<INDEL>){
        chomp $line;

        if ($line =~ /#/){
            print $line."\n";
            next;
        }
    my @F = split (/\t/, $line); 
    print join ("\t", $F[0], $F[1], $$namedata{$F[2]}, $$namedata{$F[3]}, $F[4], $F[5], $F[6], $F[7], $F[8], $F[9], $F[10], $F[11], $F[12])."\n";

    }
    close INDEL;
}
elsif ($smapfile){
  #  PLACEHOLDER
}
elsif ($bedfile){

    open (BED, $bedfile);
    foreach my $line (<BED>){
        chomp $line;

        if ($line =~ /#/){
            print $line."\n";
            next;
        }
        my @F   =  split (/\t/, $line);
        my $featcount = @F-1;
        die "Something wrong, missing ID to key for $F[0]\n" if (!$$namedata{$F[0]});

        print join ("\t", $$namedata{$F[0]}, @F[1..$featcount])."\n";

    }
    close BED

}
elsif ($agp) {
    &AGPrev($agp, $namedata);
}
else {
    die "missing parameter\n";
}

#=====================#
#  SUBROUTINES        #
#=====================#

# Read standard files
sub readfiles {
    my ($file, $type ) = @_;

    open (FILE, $file) || die "can't open file\n";

    my %data;
    foreach (<FILE>){
        chomp $_;
        next if ($_ =~ /#/);
        my @entries = split (/\t/, $_);

        if ($type eq "cmap"){
            $data{$entries[1]} = $entries[0];
        }
        elsif ($type =~ /xmap|idx/){
            $data{$entries[0]} = $entries[1];
        }
        elsif ($type =~ /agp/){
            $data{$entries[1]} = $entries[0];
        }
        else {
            die "Missing type at readfiles\n";
        }
    }
    close FILE;
    return \%data;
}


# Add Bionano idx to header whilst reading file
sub AGPrev {

    my ($file, $nkey )   = @_;

    open (FILE, $file) || die "can't open file\n";

    my %data;
    foreach (<FILE>){
        chomp $_;

        if ($_ =~ /#/){
            print "$_\n";
            next;
        }

        my @F = split (/\t/, $_);

        if ($F[4] eq "N"){
            print "$_\n";
        }
        else {
            my $idx_num = $$nkey{$F[5]};
            my $newID = $F[5]. ":IDX.$idx_num";

            $newID =~ s/(scaffold|chromosome)\:CURRENT\://;
            print join ("\t", @F[0..4], $newID, @F[6..8])."\n";
        }
    }
    close(FILE);
}
