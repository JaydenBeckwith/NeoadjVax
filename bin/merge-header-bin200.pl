#!/usr/bin/perl

use strict;

my $sample=$ARGV[0];

my $prefix = ".";

my $infile = $prefix . "/" . $sample . "/" . $sample . "_bin200_noheader.seqz";
my $outfile = $prefix . "/" . $sample . "/" . $sample . "_bin200.seqz.gz";

my $command1 = "cat seqz.header $infile | gzip > $outfile";

print STDERR "Executing $command1\n";
system("$command1");
print STDERR "Finished $command1\n";


