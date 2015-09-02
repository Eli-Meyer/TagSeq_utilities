#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if (!$ARGV[0] || $ARGV[0] eq "-h") 
	{
	print "\nCombines expression data (counts) produced by mapping multiple\n";
	print "samples against a shared database.\n";
	print "Usage: $scriptname file1 file2 .. fileN >output.tab\n"; 
	print "\toutput.tab: \tName of output file (rows=genes, columns=samples)\n\n";
	exit;
	}

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
