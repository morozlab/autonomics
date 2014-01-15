#!/usr/bin/perl -w
use strict;
use POSIX;

my ($pname) = @ARGV;
if ((not defined $pname)) {
    print "\nUsage: $0 <project_name>\n";
    print " You need to be in pipeline dir\n";
    exit 0;
}

$pname =~ s/\///;


if (-d $pname) {
  print "===================================================================\n";
  print "\t$pname\n";
  print "===================================================================\n";

  my $cmd = 'astats ' . $pname . '/' . $pname;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

  my $nr_suffix = "_blast_nr.txt";
  my $sw_suffix = "_blast_swissprot.txt";
  my $fa_suffix = "_project.fasta";
  my $fq_suffix = ".fastq";
  my $go_suffix = "_GO.txt";
  my $gocat_suffix = "_gocats.txt";
  my $pfam_suffix = "_pfam.txt";
  my $kegg_suffix = "_KEGG.txt";
  my $q_suffix = "_quantification.txt";

  my $hits = 'Sequences producing significant alignments';
  my $all = 'Query=';
  my $misses ='No hits found';

  my $prefix = " " . $pname . "/" . $pname;
  my $wc = " | wc -l";
 
  my $fasta = "./" . $pname . "/" . $pname . $fa_suffix;
  my $res = 0;
  if (-e $fasta) {
    $cmd = "grep \'>\' " . $prefix . $fa_suffix . $wc;
#    print "$cmd\n";
    $res = `$cmd`;
    chomp $res;
  }
  my $num_contigs = $res;
#  print "NUM_CONTIGS:\t$res\n";
  print "\n";
  my $nrf = "./" . $pname . "/" . $pname . $nr_suffix;
  my $nr_hits = 0;
  my $nr_misses = 0;
  $res = 0;
  if (-e $nrf) {
    $cmd = "grep \'" . $hits . "\'" . $prefix . $nr_suffix . $wc;
    $nr_hits = `$cmd`;
    chomp $nr_hits;
    $cmd = "grep \'" . $misses . "\'" . $prefix . $nr_suffix . $wc;
    $nr_misses = `$cmd`;
    chomp $nr_misses;
    $cmd = "grep " . $all . $prefix . $nr_suffix . $wc;
    $res = `$cmd`;
    chomp $res;
  }
  my $nr_total = $nr_hits + $nr_misses;

  my $p1 = 0;
  my $p2 = 0;
  if ($res) {
    $p1 = int(($nr_hits/$res)*100);
    $p2 = int(($nr_misses/$res)*100);
  }
#  print "NR: SEQS_TRIED\t$res\tSEQS_W_HITS\t$nr_hits ($p1%)\tSEQS_WO_HITS\t$nr_misses ($p2%)\tTOTAL\t$nr_total\n";
  print "BLAST_NR: SEQS_TRIED\t$res\tSEQS_W_HITS\t$nr_hits ($p1%)\tSEQS_WO_HITS\t$nr_misses ($p2%)\n";

  my $swf = "./" . $pname . "/" . $pname . $sw_suffix;
  my $sw_hits = 0;
  my $sw_misses = 0;
  $res = 0;
  if (-e $swf) {
    $cmd = "grep \'" . $hits . "\'" . $prefix . $sw_suffix . $wc;
    $sw_hits = `$cmd`;
    chomp $sw_hits;
    $cmd = "grep \'" . $misses . "\'" . $prefix . $sw_suffix . $wc;
    $sw_misses = `$cmd`;
    chomp $sw_misses;
    $cmd = "grep " . $all . $prefix . $sw_suffix . $wc;
    $res = `$cmd`;
    chomp $res;
  }
  my $sw_total = $sw_hits + $sw_misses;

  $p1 = 0;
  $p2 = 0;
  if ($res) {
    $p1 = int(($sw_hits/$res)*100);
    $p2 = int(($sw_misses/$res)*100);
  }
  print "BLAST_SWISSPROT: SEQS_TRIED\t$res\tSEQS_W_HITS\t$sw_hits ($p1%)\tSEQS_WO_HITS\t$sw_misses ($p2%)\n";
#  print "SP: SEQS_TRIED\t$res\tSEQS_W_HITS\t$sw_hits ($p1%)\tSEQS_WO_HITS\t$sw_misses ($p2%)\tTOTAL\t$sw_total\n";
  my $f = "./" . $pname . "/" . $pname . $go_suffix;
  my $num = 0;
  if (-e $f) {
      $cmd = "wc -l $f";
      $num = `$cmd`;
      my( $n, @junk ) = split( /\s+/, $num );
      $num = $n
  }
  print "GO hits: $num\n";

  $f = "./" . $pname . "/" . $pname . $gocat_suffix;
  $num= 0;
  if (-e $f) {
      $cmd = "wc -l $f";
      $num = `$cmd`;
      my( $n, @junk ) = split( /\s+/, $num );
      $num = $n
  }
  print "gocats num categories: $num\n";

  $f = "./" . $pname . "/" . $pname . $kegg_suffix;
  $num= 0;
  if (-e $f) {
      $cmd = "wc -l $f";
      $num = `$cmd`;
      my( $n, @junk ) = split( /\s+/, $num );
      $num = $n
  }
  print "KEGG hits: $num\n";

  $f = "./" . $pname . "/" . $pname . $pfam_suffix;
  my $num_hits = 0;
  my $num_files = 0;
  my $num_aa = 6 * $num_contigs;
  my $x = ceil($num_aa / 150);
  my $numf = ceil($num_aa / $x);

  if (-e $f) {
      $cmd = 'grep -v "^#" ' . $f . ' | grep -v "^$" | wc -l';
      $num_hits = `$cmd`;
      my( $n, @junk ) = split( /\s+/, $num_hits );
      $num_hits = $n;
      $cmd = "grep Copyright $f | wc -l";
      $num_files = `$cmd`;
      ( $n, @junk ) = split( /\s+/, $num_files );
      $num_files = $n;
  }
  print "pfam hits: $num_hits\n";
  print "pfam files: $num_files out of $numf\n";  

  $f = "./" . $pname . "/" . $pname . $q_suffix;
  $num= 0;
  if (-e $f) {
      $cmd = "wc -l $f";
      $res= `$cmd`;
      my( $n, @junk ) = split( /\s+/, $res );
      $num = $n;
  }
  print "quantification lines: $num\n";

=stop
  $f = "./" . $pname . "/" . $pname . $fq_suffix;
  $num= 0;
  if (-e $f) {
      $cmd = "wc -l $f";
      my $res = `$cmd`;
      my( $n, @junk ) = split( /\s+/, $res );
      $num = $n;
  }
  $num = $num/4;
  print "NUM READS: $num\n";
=cut

} else {
  print "You need to be in pipeline dir\n";
  exit 0;
}
