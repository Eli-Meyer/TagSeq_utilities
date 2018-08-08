#!/usr/bin/env perl
# written by E Meyer, elimeyer76@gmail.com
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;

Identifies PCR duplicates in a set of FASTQ sequences based on sequences of
an insert and an internal barcode in each read. Chooses the highest quality 
representative from each group for output to a non-redundant set of reads. 
Usage:   $scriptname -i input -o output <OPTIONS>
Required arguments:
	-i input		raw input reads in FASTQ format
	-o output		name for ourput file of non-redundant reads in FASTQ format
Options:
	-s barcode_start	beginning of barcode position in each read (default=1)
	-e barcode_end		end of barcode position in each read (default=4)
	-j insert_start		beginning of insert position in each read (default=10)
	-k insert_end		end of insert position in each read (default=read length)
USAGE

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SeqIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SeqIO;

# get variables from input
getopts('i:o:h:s:e:j:k');	
if (!$opt_i || !$opt_o || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
my $fastqfile = $opt_i;					# input
my $outfqfile = $opt_o;					# output
my $bstart = 1; if ($opt_s) {$bstart = $opt_s;}		# barcode start
my $bend = 4; if ($opt_e) {$bend = $opt_e;}		# barcode end
my $istart = 10; if ($opt_j) {$istart = $opt_j;}	# insert start
my $iend = "max"; if ($opt_k) {$iend = $opt_k;}		# insert end

my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");

# read input file, identify duplicates, choose highest-quality member of each group
print "Processing $fastqfile with $scriptname.\n";
print "Reading $fastqfile and identifying duplicates...\n";
system("date");
my %sh; my $scount = $initcount = $updatecount = 0;
while ($seq = $inseqs->next_seq) 
	{
	$scount++;
	$qo = new Bio::Seq::Quality(-accession_number=>$seq->display_id, -qual=>$seq->qual, -verbose=>-1);
	$qot = $qo->qual_text;
	@qoa = split(" ", $qot);
	$qid = $qo->accession_number;
	$qsum = 0; $nq = @qoa;
	foreach $q (@qoa) {$qsum += $q;}
	$meanqi = $qsum/$nq;		# get average quality score for this read
	$seqi = $seq->seq;
	if ($iend eq "max") {$iend = $seq->length;}
	$bseqi = substr($seqi, $bstart-1, ($bend-$bstart)+1);
	$iseqi = substr($seqi, $istart-1, ($iend-$istart)+1);
	$catseq = $bseqi.$iseqi;
	$sid = $seq->display_id;
	if (!exists($sh{$catseq}))	# initialize record for each unique sequence
		{
		$sh{$catseq}{"id"} = $sid;
		$sh{$catseq}{"q"} = $meanqi;
		$initcount++;
		next;
		}
	elsif ($sh{$catseq}{"q"} < $meanqi)	# update record if a higher quality copy is found
		{
		$sh{$catseq}{"id"} = $sid;
		$sh{$catseq}{"q"} = $meanqi;
		$updatecount++;
		next;
		}		
	}
# build a list of chosen reads
my %crh; my $uniquecount = 0;
foreach $r (sort(keys(%sh)))
	{
	$idi = $sh{$r}{"id"};
	$crh{$idi}++;
	$uniquecount++;
	}
print "Done.\n";
system("date");

# read input file, print chosen reads to output
print "Writing unique sequences to $outfqfile ...\n";
my $inseqs = new Bio::SeqIO(-file=>$fastqfile, -format=>"fastq");
my $outseqs = new Bio::SeqIO(-file=>">$outfqfile", -format=>"fastq");
$ocount = 0;
while ($seq = $inseqs->next_seq) 
	{
	$idi = $seq->display_id;
	if (exists($crh{$idi})) 
		{
		$outseqs->write_seq($seq); 
		$ocount++;
		}
	}
print $scount, " sequences in $fastqfile.\n";
print $ocount, " unique sequences written to $outfqfile.\n";
$pred = int(($scount-$ocount)/($scount)*1000+0.5)/10;
print $pred, "% of reads removed by this filter.\n";
system("date");
