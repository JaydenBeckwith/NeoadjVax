#!/usr/bin/perl

use strict;

my $sample=$ARGV[0];
my $sex = $ARGV[1];

my @chrs=(1..22,"X","Y");

my $conc="";
my $prefix = ".";

for(my $i=0;$i<@chrs;$i++)
{
	next if (($chrs[$i] eq 'Y') & ($sex eq 'female'));

	$conc.=$prefix . "/" . $sample . "/" . $sample . "\.chr" . $chrs[$i] . "_bin200.seqz.gz ";

	if($i==0)
	{
		system("zcat " . $prefix . "/" . $sample . "/" . $sample . "\.chr" . $chrs[$i] . "_bin200.seqz.gz | head -1 > seqz.header");
	}
}

my $outfile = $prefix . "/" . $sample . "/" . $sample . "_bin200_noheader.seqz";

my $command1 = "zcat $conc | grep -v chromosome > $outfile";

print STDERR "Executing $command1\n";
system("$command1");
print STDERR "Finished $command1\n";


