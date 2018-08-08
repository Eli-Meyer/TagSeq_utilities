#!/usr/bin/env perl
# written by E Meyer, elimeyer76@gmail.com
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Removes non-template sequences introduced at the 5' end of cDNA tags during preparation
of Tag-Seq libraries. These regions are identified as the 3'-most occurence of the GGG
motif within a user-defined window (based on RNA oligo sequences used in library preparation).
Reads lacking a GGG motif in the chosen region are discarded.
Usage:   $scriptname -i input -o output <OPTIONS>
Required arguments:
	-i input	raw input reads in FASTQ format
	-o output	name for ourput file of HQ reads in FASTQ format
Options:
	-b beginning	defines the beginning of the range in which to search for GGG. (default=1)
	-e end		defines the end of range in which to search for GGG. (default=8)

USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:o:b:e:h');	# in this example a is required, b is optional, h is help
if (!$opt_i ||!$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
my $seqfile = $opt_i;
my $outfile = $opt_o;
my $rstart = 1; if ($opt_s) {$rstart = $opt_s;}	# beginning of region to search
my $rend = 8; if ($opt_e) {$rend = $opt_e;}	# end of region to search

# loop through fastq file and record presence and position of motif in each read
print "Processing $seqfile with $scriptname.\n";
print "Reading $seqfile and searching for GGG motif in the region $rstart - $rend...\n";
system("date");
open (IN, $seqfile);
my $line = 0;
while(<IN>)
	{
	chomp;
	$line++; 
	if ($line eq 1) 
		{
		$sid = $_;
		$incount++;
		next;
		}
	elsif ($line eq 2)
		{
		$ssi = $_;
		$regseq = substr($ssi, $rstart-1, ($rend-$rstart+1));
		$found = 0;
		for ($a=$rend-3; $a>=$rstart; $a--)
			{
			$subseq = substr($regseq, $a, 3);
			if ($subseq eq "GGG")
				{
				$cuti = $a+3+1;
				$found++;
				last;
				}
			}
		if ($found eq 0)
			{
			$nomotif++;
			next;
			}
		$motif++;
		$sh{$sid} = $cuti;
#		print $sid, "\t", $regseq, "\t", $cuti, "\t";
#		print $ssi, "\t", substr($ssi, $cuti-1, (length($ssi)-$cuti)), "\n";
		next;
		}
	elsif ($line eq 4) {$line = 0;}
	}
close(IN);
print "Done.\n";
system("date");

# read in fastq file again and print trimmed reads to output
print "Printing trimmed reads to $outfile ...\n";
open (IN, $seqfile);
open (OUT, ">$outfile");
my $line = 0;
while(<IN>)
	{
	chomp;
	$line++; 
	if ($line eq 1) 
		{
		$sid = $_;
		}
	if (exists($sh{$sid}))
		{
		if ($line eq 1 || $line eq 3)
			{
			print OUT $_, "\n";
			}
		if ($line eq 2 || $line eq 4)
			{
			$ssi = $_;
			$leni = length($ssi);
			$subseq = substr($ssi, ($sh{$sid}-1), ($leni-$sh{$sid}));
			print OUT $subseq, "\n"; 
			}
		if ($line eq 1) {$outcount++;}
		}
	if ($line eq 4) {$line = 0;}
	}
close(IN);
close(OUT);
print "Done.\n";
system("date");
print $incount, " sequences in $seqfile.\n";
print $motif, " of these had a GGG motif in the region $rstart - $rend.\n";
$pred = int($nomotif/$incount*1000+0.5)/10;
print $nomotif, " of these lacked this motif in this region ($pred %).\n";
print $outcount, " sequences written to $outfile.\n";
