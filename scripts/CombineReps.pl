#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Combines counts data from RNASeq libraries prepared from the same biological
sample (e.g. replicate libraries or sequencing runs).
Usage: $scriptname input_1 input_2 ... input_n > output
Where:
	input_1:	first counts file (tab delimited text; rows=genes, columns=samples)
	...
	input_n:	final counts file to be combined
	output:		a name for the output file (same format as input)
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- loop through input files and record counts for each sample and gene
foreach $argi (0..$#ARGV)
	{
	open(IN,$ARGV[$argi]);
	$rowcount=0;
	while(<IN>)
		{
		chomp;
		$rowcount++; 
		@cols = split("\t", $_);
		$sumh{$cols[0]}{$argi}+=$cols[1];
		}
	}

# -- loop through expressed genes and report the sum of all counts for each gene
print "\t";
foreach $argi (0..$#ARGV) {print $argi, "\t";}
print "\n";

@cgenes = sort(keys(%sumh));
foreach $gene (@cgenes)
	{
	print $gene, "\t";
	foreach $argi (0..$#ARGV) {print $sumh{$gene}{$argi}, "\t";}
	print "\n";
	}
