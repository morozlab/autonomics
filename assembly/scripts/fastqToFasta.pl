use warnings;
use strict;
use Bio::SeqIO;

my $usage = "fastqToFasta.pl <fastqFile> <desired prefix (including path)>\n";

my $file = shift or die $usage;
my $prefix = shift or die $usage;

my $seqio = new Bio::SeqIO(-format=>"fastq", -file=>"$file");

open(SEQFILE, ">$prefix\.fasta");
open(QUALFILE, ">$prefix\.fasta.qual");

while(my $seq = $seqio->next_seq){
	my $id = $seq->display_id;
	if ( defined( $seq->desc ) ) {
		$id .= " " . $seq->desc;
	}
	$id =~ s/(^\s+|\s+$)//g;
	print SEQFILE ">$id\n";
	print QUALFILE ">$id\n";
	print SEQFILE $seq->seq . "\n";
	print QUALFILE $seq->qual_text . "\n";
}

close(SEQFILE);
close(QUALFILE);


