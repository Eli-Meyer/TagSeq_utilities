#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Filters a FASTQ file to remove sequencings containing homopolymer repeats (HR) 
(e.g. polyA tails) longer than the specified threshold. 
Usage:	$scriptname -i input -n max_length -o output
Required arguments:
	-i input	a file of short reads to be filtered (FASTQ)
	-n max_length	maximum length of homopolymer repeat allowed in output
	-o output	a name for the output file of sequences passing this filter (FASTQ)
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:n:o:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_n || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

my $seqfile = $opt_i;		# raw reads, fastq format
my $critlen = $opt_n;		# critical length
my $outfile = $opt_o;		# name for output file, fastq format
my @vals = qw{A C G T N};
my @stra = ();
foreach $v (@vals)
	{
	$ssi = $v x $critlen;
	push (@stra, $ssi);
	}

# loop through fastq file and print out only those passing filter
print "Processing $seqfile with $scriptname ...\n";
system("date");
open (IN, $seqfile);
open (OUT, ">$outfile");
my $switch = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) 
		{
		$thisname = $_;
		$nextup=1;
		$switch=0;
		}
	else
		{
		if ($nextup==1)
			{
			$switch=1;
			foreach $s (@stra) 
				{
				if ($_ =~ /$s/) {$switch=0;}
				}			
			$nextup=0;
			if ($switch==0) {$fail++;}
			elsif ($switch==1) 
				{
				print OUT $thisname, "\n";
				$pass++;
				}
			}
		if ($switch>0) {print OUT $_, "\n";}
		}
	}
close(IN);
print $fail, " reads failed\n";
print $pass, " reads passed\n";
system("date");
