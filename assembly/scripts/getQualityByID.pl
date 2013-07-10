use warnings;
use strict;
use Bio::SeqIO;

my $usage =
"getQualityByID.pl <fasta file> <file of sequence IDs, one per line> <file format>\n";

my $qualFile = shift or die $usage;
my $list     = shift or die $usage;
my $format   = shift or die $usage;

my %idHash = ();
open( LIST, "<$list" );

while (<LIST>) {
	chomp();
	$_ =~ s/(^\s+|\s+$)//g;
	my @array = split( /\s+/, $_ );
	$idHash{ $array[0] } = 1;
}
if ( $format eq "fasta" ) {

	#open the quality file
	open( QUAL, "<$qualFile" ) or die "Could not open quality file\n";
	my $qualID  = "";
	my $quality = "";
	while (<QUAL>) {
		chomp();
		my $line = $_;
		if ( $line =~ /^>/ ) {
			if ( $qualID ne "" ) {

				#check if we should print this quality score
				if ( exists( $idHash{$qualID} ) ) {
					print ">$qualID\n$quality\n";
				}
				$quality = "";
			}
			my @array = split( /\s+/, $line );
			$line   = $array[0];
			$qualID = $line;
			$qualID =~ s/^>//;
			$qualID =~ s/(^\s+|\s+$)//g;
		}
		else {
			if ( $quality eq "" ) {
				$quality = $line . " ";
			}
			else {
				$quality .= $line . " ";
			}
		}
	}
	if ( exists( $idHash{$qualID} ) ) {
		print ">$qualID\n$quality\n";
	}
}
elsif ( $format eq "fastq" ) {
	my $seqio = new Bio::SeqIO( -format => "fastq", file => "$qualFile" );
	while ( my $seq = $seqio->next_seq ) {
		my $id = $seq->display_id;
		if ( defined( $seq->desc ) ) {
			$id .= " " . $seq->desc;
		}
		$id =~ s/(^\s+|\s+$)//g;
		if ( exists( $idHash{$id} ) ) {
			print ">" . $seq->display_id . "\n";
			print $seq->qual_text . "\n";
		}
	}
}
