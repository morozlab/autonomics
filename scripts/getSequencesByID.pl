use warnings;
use strict;
use Bio::SeqIO;

my $usage =
"getSequencesByID.pl <sequence file> <file of sequence IDs, one per line> <file format> <min-length> <type of file>\n";

my $file   = shift or die $usage;
my $list   = shift or die $usage;
my $format = shift or die $usage;
my $minLength = shift;
my $type   = shift;

if(!(defined($minLength))){
	$minLength = 0;
}

my %idHash = ();
open( LIST, "<$list" );

my $length = 0;
while (<LIST>) {
	chomp();
	$_ =~ s/(^\s+|\s+$)//g;
	$idHash{$_} = 1;
}

my $seqio = new Bio::SeqIO( -format => "$format", -file => "$file" );

while ( my $seq = $seqio->next_seq ) {
	my $id = $seq->display_id;
	#if ( defined( $seq->desc ) ) {
	#	$id .= " " . $seq->desc;
	#}
	$id =~ s/(^\s+|\s+$)//g;
	if ( defined($type) ) {
		if ( $type eq "velvet_genome" ) {
			my @array = split( /\|/, $id );
			$id = $array[0];
		}
	}

	if ( exists( $idHash{$id} ) ) {
		if(length($seq->seq) > $minLength){
		    print ">" . $seq->display_id . "\n";
		    print $seq->seq . "\n";
		}
	}
}

