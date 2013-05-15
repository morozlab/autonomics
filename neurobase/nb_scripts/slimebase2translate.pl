use Bio::Seq;
use Bio::SeqUtils;
use strict;

my $usage = "usage: translateFasta.pl <single sequence passed via amfphp gateway>";
my $seq = shift or die $usage . "\n";
my $seqobj = new Bio::Seq(-display_id=>"1", -seq=>$seq);
my $seqline = $seqobj->seq;
$seqline =~ tr/(X|\#)/N/;
$seqobj->seq($seqline);
my @proteinSeqs = Bio::SeqUtils->translate_6frames($seqobj);
foreach my $protein (@proteinSeqs){
   print $protein->seq . "\n";
}


