#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Analyzes a set of alignments, filtering out weak or uninformative alignments
and counting the number of reads matching each gene in the reference. 

Matches may be counted as:
(1) the number of reads matching each sequence in the reference
(2) the number of reads matching each gene, from a user defined list
(3) the number of reads matching each component or subcomponent in a Trinity assembly

NOTE: make sure that when a read matches multiple reference sequences (ambigous)
your mapper reports all of these, or at least all alignments as strong as the best alignment. 
e.g. with SHRiMP you could use the flag --strata. This is NOT the default behavior for some 
mappers, but is required to exclude ambiguous matches before further analysis.

Usage:  $scriptname -i input -m matches -o output <options>
Required arguments:
	-i input	Output from any short read mapper, in SAM format.
	-m matches	Minimum number of matching bases required to consider an alignment valid. 
	-o output	A name for the filtered output (SAM format). 
Options:
	-p option	1: Report the number of reads matching each reference sequence
			in a separate output files "counts.tab". 
			0: Don't produce this file (default).
	-r method	The method used for counting matches.
			s: (default) count the number of reads uniquely assigned to a single reference sequence.
			t: count the number of reads uniquely assigned to each gene, as defined by the 
			   component-subcomponent-transcript structure in denovo transcriptome assemblies
			   produced by the Trinity assembler. This requires that your reference sequences are 
			   named in the style of Trinity assemblies, e.g. comp0_c0_seq1. (see option -t)
			g: count the number of reads uniquely assigned to each gene, as defined by a 
			   user-supplied gene list (see option -g)
	-t option	Options for method "-r t". Choose whether to count genes at the level of
			components ("-t c")(the default) or sub-components ("-t s").
	-g gene_list	Required for method "-r g". The name of a gene list (tab-delimited text) formatted as
			  sequence1	gene1
			  sequence2	gene1
	-l length	Minimum length of aligned region (match, mismatch, + gaps) required to consider 
			an alignment valid. Only relevant if your mapper uses local alignment. Default (for
			global alignments) is set equal to -m. 
USAGE
if ($#ARGV < 3 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('i:m:o:p:r:t:g:l:h');	# in this example a is required, b is optional, h is help
if (!$opt_i || !$opt_m || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_p) {$cprint = $opt_p;} else {$cprint = 0;}
if ($opt_t) {$countlevel = $opt_t;} else {$countlevel = "c";}
if ($opt_r) {$method = $opt_r;} else {$method = "s";}
if ($opt_r eq "g" && !$opt_g) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
if ($opt_g)	{$amblist = $opt_g;}
if ($opt_l) {$athd = $opt_l;} else {$athd = $opt_m;}
my $infile = $opt_i;
my $mthd = $opt_m;
my $outfile = $opt_o;
$ambig = $tooshort = 0;

if ($method eq "s")
{
# read in sam file output and build a hash, counting raw mappings
open(IN, $infile);
my %maph;
while(<IN>)
	{
	if ($_ =~ /^@/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$ncols = @cols;

	if ($cols[2] eq "*") {next;}
	$rawmap++;

# -- extract alignment length 
        $numstr = $cols[5];
        @chars = split("M", $numstr);
	$aligni = 0;
	foreach $c (@chars)
		{
		$c =~ s/.+\D//g;
		if ($c > 0) {$aligni+=$c;}
		}
# -- extract mismatches and apply alignment length and mismatch thresholds
	$mismatchi = 0;
	foreach $flag (@cols[11..$ncols])
		{
		if ($flag =~ /NM\:i/)
			{
			$flag =~ s/NM\:i\://;
			$mismatchi = $flag;
			}
		}
	$matchi = $aligni - $mismatchi;
	if ($matchi < $mthd) {$tooweak++; next;}
	if ($aligni < $athd) {$tooshort++; next;}

# -- add extracted data to hash
	$maph{$cols[0]}{$cols[2]}{"count"}++;
	$maph{$cols[0]}{$cols[2]}{"string"} = $_;
	$maph{$cols[0]}{$cols[2]}{"align"} = $aligni;
	$maph{$cols[0]}{$cols[2]}{"match"} = $matchi;
	}

# count number of reads with one or more matches passing thresholds
@mra = keys(%maph);
$nmr = @mra;

# select final set of unique mappings
open(OUT, ">$outfile");
open(IN, $infile);
while(<IN>)
	{
	chomp;
	if ($_ =~ /^@/) {print OUT $_, "\n";}
	else {last;}
	}
my %refch;
foreach $r (@mra)
	{
	%rh = %{$maph{$r}};
	@ma = keys(%rh);
	$nma = @ma;
	if ($nma>1)
		{
		@moa = sort{$rh{$b}{"match"}<=>$rh{$a}{"match"}}(keys(%rh));
		if ($rh{$moa[0]}{"match"}==$rh{$moa[1]}{"match"}) 
			{
			$ambig++; 
			next;
			}
		}
	if ($nma == 1) {@moa = @ma;}
	$unimap++;
	$refch{$moa[0]}++;
	print OUT $rh{$moa[0]}{"string"}, "\n";
	}

if ($cprint eq 1)
{
open(CTS, ">counts.tab");
@sref = sort{$refch{$b}<=>$refch{$a}}(keys(%refch));
foreach $s (@sref)
	{
	print CTS $s, "\t", $refch{$s}, "\n";
	}
}

print $rawmap, " raw mappings altogether.\n";
print $nmr, " reads had one or more matches\n";
print $tooshort, " excluded for short alignments\n";
print $tooweak, " excluded for weak matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $unimap, " unique mappings remained\n";
}

elsif ($method eq "g")
{
# read in list of functionally equivelent reference sequences (FERS) and build a hash
%fersh = ();
open(IN, $amblist);
while(<IN>)
	{
	chomp;
	($seq, $gene) = split("\t", $_);
	$fersh{$seq} = $gene;
	}
close(IN);

# read in sam file output and build a hash, counting raw mappings
$ambig = $tooshort = 0;
open(IN, $infile);
my %maph = ();
while(<IN>)
	{
	if ($_ =~ /^@/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$ncols = @cols;
	if ($cols[2] eq "*") {next;}
	$rawmap++;

        $numstr = $cols[5];
        @chars = split("M", $numstr);

# -- extract alignment length 
        $numstr = $cols[5];
        @chars = split("M", $numstr);
	$aligni = 0;
	foreach $c (@chars)
		{
		$c =~ s/.+\D//g;
		if ($c > 0) {$aligni+=$c;}
		}
# -- extract mismatches and apply alignment length and mismatch thresholds
	$mismatchi = 0;
	foreach $flag (@cols[11..$ncols])
		{
		if ($flag =~ /NM\:i/)
			{
			$flag =~ s/NM\:i\://;
			$mismatchi = $flag;
			}
		}
	$matchi = $aligni - $mismatchi;
	if ($matchi < $mthd) {$tooweak++; next;}
	if ($aligni < $athd) {$tooshort++; next;}

# -- check the FERS list for each sequence to decide how to count the reads
	$assign = "";
	if(exists($fersh{$cols[2]}))
		{
		$assign = $fersh{$cols[2]};
		}
	else	{$assign = $cols[2];}
# -- check if a mapping already exists for this read-ref pair
# -- and keep the best mapping if so
	if(exists($maph{$cols[0]}{$assign}))
		{
		$prevmatch = $maph{$cols[0]}{$assign}{"match"};
		$thismatch = $matchi;
		if ($prevmatch >= $thismatch) {next;}
		}
# -- add extracted data to hash
	$maph{$cols[0]}{$assign}{"count"}++;
	$maph{$cols[0]}{$assign}{"string"} = $_;
	$maph{$cols[0]}{$assign}{"align"} = $aligni;
	$maph{$cols[0]}{$assign}{"match"} = $matchi;
	}

# count number of reads with one or more matches passing thresholds
@mra = keys(%maph);
$nmr = @mra;

# select final set of unique mappings
open(OUT, ">$outfile");
open(IN, $infile);
while(<IN>)
	{
	chomp;
	if ($_ =~ /^@/) {print OUT $_, "\n";}
	else {last;}
	}
my %refch;
foreach $r (@mra)
	{
	%rh = %{$maph{$r}};
	@ma = keys(%rh);
	$nma = @ma;
	if ($nma>1)
		{
		@moa = sort{$rh{$b}{"match"}<=>$rh{$a}{"match"}}(keys(%rh));
		if ($rh{$moa[0]}{"match"}==$rh{$moa[1]}{"match"} && $moa[0] ne $moa[1]) 
			{
			$ambig++; 
			next;
			}
		}
	if ($nma == 1) {@moa = @ma;}
	$unimap++;
	$refch{$moa[0]}++;
	print OUT $rh{$moa[0]}{"string"}, "\n";
	}

if ($cprint eq 1)
{
open(CTS, ">counts.tab");
@sref = sort{$refch{$b}<=>$refch{$a}}(keys(%refch));
foreach $s (@sref)
	{
	print CTS $s, "\t", $refch{$s}, "\n";
	}
}
print $rawmap, " raw mappings altogether.\n";
print $nmr, " reads had one or more matches\n";
print $tooshort, " excluded for short alignments\n";
print $tooweak, " excluded for weak matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $unimap, " unique mappings remained\n";
}

elsif ($method eq "t")
{
# read in sam file output and build a hash, counting raw mappings
open(IN, $infile);
my %maph;
while(<IN>)
	{
	if ($_ =~ /^@/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$ncols = @cols;

	if ($cols[2] eq "*") {next;}
	$rawmap++;

        $numstr = $cols[5];
        @chars = split("M", $numstr);

# -- extract alignment length 
        $numstr = $cols[5];
        @chars = split("M", $numstr);
	$aligni = 0;
	foreach $c (@chars)
		{
		$c =~ s/.+\D//g;
		if ($c > 0) {$aligni+=$c;}
		}
# -- extract mismatches and apply alignment length and mismatch thresholds
	$mismatchi = 0;
	foreach $flag (@cols[11..$ncols])
		{
		if ($flag =~ /NM\:i/)
			{
			$flag =~ s/NM\:i\://;
			$mismatchi = $flag;
			}
		}
	$matchi = $aligni - $mismatchi;
	if ($matchi < $mthd) {$tooweak++; next;}
	if ($aligni < $athd) {$tooshort++; next;}

# -- apply the count level decision
	$refname = $cols[2];
	if ($countlevel eq "c")		{$refname =~ s/_.+//g;}
	elsif ($countlevel eq "s")	{$refname =~ s/_seq.+//g;} 
# -- add extracted data to hash
	$maph{$cols[0]}{$refname}{"count"}++;
	$maph{$cols[0]}{$refname}{"string"} = $_;
	$maph{$cols[0]}{$refname}{"align"} = $aligni;
	$maph{$cols[0]}{$refname}{"match"} = $matchi;
	}

# count number of reads with one or more matches passing thresholds
@mra = keys(%maph);
$nmr = @mra;

# select final set of unique mappings
open(OUT, ">$outfile");
open(IN, $infile);
while(<IN>)
	{
	chomp;
	if ($_ =~ /^@/) {print OUT $_, "\n";}
	else {last;}
	}
my %refch;
foreach $r (@mra)
	{
	%rh = %{$maph{$r}};
	@ma = keys(%rh);
	$nma = @ma; $status = 0;
	if ($nma>1)
		{
		@moa = sort{$rh{$b}{"match"}<=>$rh{$a}{"match"}}(keys(%rh));
		for ($m=0; $m<@moa; $m++)
			{
			if ($status > 0) {next;}	
			if ($rh{$moa[$m]}{"match"}==$rh{$moa[$m+1]}{"match"}) 
				{
				if ($moa[$m] ne $moa[$m+1])
					{
					$ambig++;
					$status++;
					}
				} 
			}
		}
	if ($nma == 1) {@moa = @ma;}
	if ($status > 0) {next;}
	$unimap++;
	$refch{$moa[0]}++;
	print OUT $rh{$moa[0]}{"string"}, "\n";
	}

if ($cprint eq 1)
{
open(CTS, ">counts.tab");
@sref = sort{$refch{$b}<=>$refch{$a}}(keys(%refch));
foreach $s (@sref)
	{
	print CTS $s, "\t", $refch{$s}, "\n";
	}
}

print $rawmap, " raw mappings altogether.\n";
print $nmr, " reads had one or more matches\n";
print $tooshort, " excluded for short alignments\n";
print $tooweak, " excluded for weak matches\n";
print $ambig, " excluded for ambiguous matches\n";
print $unimap, " unique mappings remained\n";
}
