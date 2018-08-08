#!/usr/bin/env perl
# written by E Meyer, elimeyer76@gmail.com
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Combines expression data (counts) produced by mapping multiple
samples against the same set of reference sequences.
Usage: $scriptname file_1 file_2 ... file_n > output_file
Where:
	files 1-n:	Gene expression data (counts output from SAMFilterByGene.pl) for each sample
	output_file:	a name for the output; tab-delimited text with genes as rows and samples as columns
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
