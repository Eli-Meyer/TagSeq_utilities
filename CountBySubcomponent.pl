#!/usr/bin/perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

$scriptname=$0; $scriptname =~ s/.+\///g;
if ($#ARGV != 3 || $ARGV[0] eq "-h") 
	{
	print "\nFilters alignments produced by SHRiMP's gmapper and probcalc to exclude weak or short,\n";
	print "partial matches based on user specified criteria.\n";
	print "If your reference is a Trinity de novo assembly, this script may be appropriate.\n";
	print "Ambiguous mapping is evaluated at the level of subcomponents: if a read matches two sequences\n";
	print "from the same subcomponent, that read is counted. If a read matches two sequences from different\n";
	print "subcomponents, that read is discarded.\n\n";
	print "Usage: $scriptname input_alignments min_matching min_align output_counts\n"; 
	print "Where: input_alignments:\tinput file, the output from gmapper -> probcalc.\n";
	print "\tmin_matching:\t\tMinimum number of matching bases required in each alignment.\n";
	print "\tmin_align:\t\tMinimum length of each alignment required.\n";
	print "\toutput_counts:\t\tA name for the output file documenting coverage for each reference sequence.\n\n";
	exit;
	}

my $infile = $ARGV[0];		# output from probcalc
my $mthd = $ARGV[1];		# min number of matching bases
my $athd = $ARGV[2];		# min number of aligned bases
my $outfile = $ARGV[3];		# name for output

# -- read in gmapper output and record passing matches by read 
$ambig = 0;
print "Recording passing matches for each read...\n";
open(IN, $infile);
while(<IN>)
	{
	if ($_ =~ /^#/) {next;}
	chomp;
	$rowi = $_;
	$rowi =~ s/^>//;
	@cols = split("\t", $rowi);
	$rawmap++;
	$al = $cols[6]-$cols[5]+1;
        $numstr = $cols[9];
        @chars = split("", $numstr);
        my $mmno = 0;
        foreach $c (@chars)
                {
		if ($insw == 0)	{
		if ($c =~ /\d/) {next;}
		if ($c =~ /[xACGT-]/i) {$mmno++; next;}
		if ($c =~ /\(/) {$insw = 1; next;}
				}
		if ($insw == 1)	{
		if ($c =~ /\d/) {$mmno++; next;}
		if ($c =~ /\)/) {$insw = 0; next;}
				}	
                }
	$matchi = $al - $mmno;
	if ($al < $athd || $matchi < $mthd) {$tooshort++; next;}
	$passh{$cols[0]}++;
#	print $cols[0], "\t", $cols[1], "\t", $cols[8], "\n";
	if (exists($rmh{$cols[0]}))
		{
		if ($cols[8] >= $rmh{$cols[0]}{"score"})
			{
			$rmh{$cols[0]}{"matches"} = $rmh{$cols[0]}{"matches"}." ".$cols[1];
			}
		}
	else
		{
		$rmh{$cols[0]}{"matches"} = $cols[1];
		}
	}
print "Done.\n";
# -- count number of reads with one or more matches
@pra = keys(%passh);
$pmaps = @pra;

# -- loop through all matches for each read and evaluate whether more than one subcomponent was hit
print "Checking for ambiguity...\n";
foreach $r (sort(keys(%rmh)))
	{
	$longenough++;
	%igih = ();
	@ama = split(" ", $rmh{$r}{"matches"});
	if (@ama>1) {$multcount++;}
	foreach $a (@ama)
		{
		$a =~ s/_seq\d+//;
		$igih{$a}++;
		}
	@iga = sort(keys(%igih));
#	print $r, "\t", "@ama", "\n", "@iga", "\n";

# -- exclude reads mapping to multiple subcomponents
	$niga = @iga;
	if ($niga > 1) {$ambig++; next;}

# -- count passing reads by subcomponent
	$unambig++;
	$ch{$iga[0]}++;
	if (@ama==1) {$singref++;}
	}
print "Done.\n";

# -- print out counts by subcomponent
print "Writing output...\n";
open(OUT, ">$outfile");
foreach $c (sort(keys(%ch)))
	{
	print OUT $c, "\t", $ch{$c}, "\n";
	}
close(OUT);
print "Done.\n";
print "Process completed.\n";

# -- print a summary
print $rawmap, " raw mappings in input.\n";
print $tooshort, " mappings excluded for being too short.\n";
print $longenough, " mappings were long enough to consider.\n";
print $multcount, " reads mapped to multiple transcript equally well.\n";
print $ambig, " excluded for matching different subcomponents equally well.\n";
print $singref, " showed unambiguous matches to a single transcript.\n";
print $unambig, " showed unambiguous matches to a single subcomponent.\n";
