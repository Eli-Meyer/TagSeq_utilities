-------------------------
CombineExpression.pl
-------------------------

------------------------------------------------------------
CombineExpression.pl
Combines expression data (counts) produced by mapping multiple
samples against the same set of reference sequences.
Usage: CombineExpression.pl file_1 file_2 ... file_n > output_file
Where:
	files 1-n:	Gene expression data (counts output from SAMFilter.pl) for each sample
	output_file:	a name for the output; tab-delimited text with genes as rows and samples as columns
------------------------------------------------------------

-------------------------
HRFilterFastq.pl
-------------------------

------------------------------------------------------------
HRFilterFastq.pl
Filters a FASTQ file to remove sequencings containing homopolymer repeats (HR) 
(e.g. polyA tails) longer than the specified threshold. 
Usage:	HRFilterFastq.pl -i input -n max_length -o output
Required arguments:
	-i input	a file of short reads to be filtered (FASTQ)
	-n max_length	maximum length of homopolymer repeat allowed in output
	-o output	a name for the output file of sequences passing this filter (FASTQ)
------------------------------------------------------------

-------------------------
QualFilterFastq.pl
-------------------------

------------------------------------------------------------
QualFilterFastq.pl

Removes reads containing too many low quality basecalls from a set of short sequences
Output:  high-quality reads in FASTQ format
Usage:   QualFilterFastq.pl -i input -m min_score -x max_LQ -o output
Required arguments:
	-i input	raw input reads in FASTQ format
	-m min_score	quality scores below this are considered low quality (LQ)
	-x max_LQ	reads with more than this many LQ bases are excluded
	-o output	name for ourput file of HQ reads in FASTQ format
------------------------------------------------------------

-------------------------
RemovePCRDups.pl
-------------------------

------------------------------------------------------------
RemovePCRDups.pl

Identifies PCR duplicates in a set of FASTQ sequences based on sequences of
an insert and an internal barcode in each read. Chooses the highest quality 
representative from each group for output to a non-redundant set of reads. 
Usage:   RemovePCRDups.pl -i input -o output <OPTIONS>
Required arguments:
	-i input		raw input reads in FASTQ format
	-o output		name for ourput file of non-redundant reads in FASTQ format
Options:
	-s barcode_start	beginning of barcode position in each read (default=1)
	-e barcode_end		end of barcode position in each read (default=4)
	-j insert_start		beginning of insert position in each read (default=10)
	-k insert_end		end of insert position in each read (default=read length)
------------------------------------------------------------

-------------------------
SAMFilterByGene.pl
-------------------------

------------------------------------------------------------
SAMFilterByGene.pl
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

Usage:  SAMFilterByGene.pl -i input -m matches -o output <options>
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
------------------------------------------------------------

-------------------------
TagTrimmer.pl
-------------------------

------------------------------------------------------------
TagTrimmer.pl

Removes non-template sequences introduced at the 5' end of cDNA tags during preparation
of Tag-Seq libraries. These regions are identified as the 3'-most occurence of the GGG
motif within a user-defined window (based on RNA oligo sequences used in library preparation).
Reads lacking a GGG motif in the chosen region are discarded.
Usage:   TagTrimmer.pl -i input -o output <OPTIONS>
Required arguments:
	-i input	raw input reads in FASTQ format
	-o output	name for ourput file of HQ reads in FASTQ format
Options:
	-b beginning	defines the beginning of the range in which to search for GGG. (default=1)
	-e end		defines the end of range in which to search for GGG. (default=8)

------------------------------------------------------------

