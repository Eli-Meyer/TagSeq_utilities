#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Combines gene expression results from multiple samples, producing an output of gene
counts (the number of reads mapping to each gene from each sample) suitable for
statistical analysis of gene expression.

Usage: $scriptname file_1 file_2 ... file_n > output_file
Where:
	files 1-n:	nucleotide frequencies (output from SAMBasecaller.pl) for each sample
	output_file:	a name for the output; tab-delimited text
USAGE
if ($#ARGV < 1 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my %bigh;

print "\t";
foreach $argi (0..$#ARGV)
	{
	print $ARGV[$argi], "\t";
	open(TAB, $ARGV[$argi]);
	while(<TAB>)
		{
		chomp;
		@cols = split("\t", $_);
		$bigh{$cols[0]}{$argi} = $cols[1];		
		}
	}
print "\n";

foreach $r (sort(keys(%bigh)))
	{
	print $r, "\t";
	foreach $argi (0..$#ARGV)
		{
		if(exists($bigh{$r}{$argi}))
			{print $bigh{$r}{$argi}, "\t";
			}
		else	{print 0, "\t";}
		}
	print "\n";
	}
