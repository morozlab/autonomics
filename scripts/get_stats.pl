#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my $ASSEMBLE = 'assemble';

my $QUANTIFICATION = 'quantification';
my $GO = 'go';
my $KEGG = 'kegg';
my $PFAM = 'pfam';
my $BLAST_NR = 'blast_nr';
my $SWISS = 'blast_swissprot';
my $PANTHER = 'panther';

my $READS_EXT =  '.fastq';
my $ASSEMBLE_EXT = '_project.fasta';
my $QUANTIFICATION_EXT = '_quantification.txt';
my $GO_EXT = '_GO.txt';
my $GOCATS_EXT = '_gocats.txt';
my $KEGG_EXT = '_KEGG.txt';
my $PFAM_EXT = '_pfam.txt';
my $BLAST_NR_EXT = '_blast_nr.txt';
my $PROTEINS_EXT = '_project.fasta_translated.fa';
my $SWISS_EXT = '_blast_swissprot.txt';
my $PANTHER_EXT = '_panther';

my ($job_type, $project_name) = @ARGV;
if ((not defined $job_type) || (not defined $project_name)) {
    print "\nUsage: $0 <job_type> <project_name>\n";
    exit 0;
}

my $cmd;
my $data_file;
my $stats_file = $project_name . ".stats";
open STATS, ">>$stats_file" or die "Unable to open $stats_file:";

if ($job_type eq $ASSEMBLE) {
  $data_file = $project_name . $ASSEMBLE_EXT;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "sdb2 $data_file";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)\t(.*)/;
      print STATS "num_assembled_seqs\t$1\n";
      print STATS "num_assembled_bases\t$2\n";
    }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $data_file = $project_name . $READS_EXT;

  my $total_bases_and_new_lines = 0;
  my $total_new_lines = 0;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = 'grep -B 1 "^+$" '. $data_file . ' | grep -v "\-\-" | grep -v "+" | wc -m';
    print "$cmd\n";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      print "num_bases_and_new_lines\t$1\n";
      $total_bases_and_new_lines = $1;
    }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $cmd = "wc -l $data_file";
  print "$cmd\n";
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*) (.*)/;
    print "num_new_lines\t$1\n";
    $total_new_lines = $1;
  }
  my $num_raw_bases = $total_bases_and_new_lines - $total_new_lines;
  print STATS "num_raw_bases\t$num_raw_bases\n";
  print STATS "num_raw_reads\t$total_new_lines\n";
} elsif ($job_type eq $BLAST_NR) {
  $data_file = $project_name . $BLAST_NR_EXT;
  my $total_num_queries = 0;
  my $num_misses = 0;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "grep Query= $data_file | wc -l";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      $total_num_queries = $1;
      print "total_num_queries = $total_num_queries\n";
      }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $cmd = 'grep "\*\*\*\*" ' . $data_file . ' | wc -l';
  print "$cmd\n";
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    $num_misses = $1;
    print "num_misses = $num_misses\n";
  }
  my $num_successful_queries = $total_num_queries - $num_misses;
  print STATS "num_transcripts_with_blast_nr_hits\t$num_successful_queries\n";
  $cmd = 'grep ">" ' . $data_file . ' | sort -u | wc -l';
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    print STATS "num_blast_nr_hits\t$1\n";
  }
} elsif ($job_type eq $SWISS) {
  $data_file = $project_name . $SWISS_EXT;
  my $total_num_queries = 0;
  my $num_misses = 0;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "grep Query= $data_file | wc -l";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      $total_num_queries = $1;
      print "total_num_queries = $total_num_queries\n";
      }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $cmd = 'grep "\*\*\*\*" ' . $data_file . ' | wc -l';
  print "$cmd\n";
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    $num_misses = $1;
    print "num_misses = $num_misses\n";
  }
  my $num_successful_queries = $total_num_queries - $num_misses;
  print STATS "num_transcripts_with_blastswissprot_hits\t$num_successful_queries\n";

  $cmd = 'grep ">" ' . $data_file . ' | sort -u | wc -l';
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    print STATS "num_blastswissprot_hits\t$1\n";
  }
} elsif ($job_type eq $GO) {
  $data_file = $project_name . $GO_EXT;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "grep GO: $data_file | wc -l";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      print STATS "num_go_hits\t$1\n";
      }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $cmd = 'cut -f1,1 -d " " ' . $data_file . ' | sort -u | wc -l';
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    print STATS "num_transcripts_with_go_hits\t$1\n";
  }
} elsif ($job_type eq $KEGG) {
  $data_file = $project_name . $KEGG_EXT;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "wc -l $data_file";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      print STATS "num_kegg_hits\t$1\n";
      }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }
  $cmd = 'cut -f1,1 -d " " ' . $data_file . ' | sort -u | wc -l';
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    print STATS "num_transcripts_with_kegg_hits\t$1\n";
  }
} elsif ($job_type eq $PFAM) {
  $data_file = $project_name . $PFAM_EXT;
  if ((-e $data_file) && (not -z $data_file)) {
    $cmd = "grep PF $data_file | wc -l";
    open (GR, "$cmd |") or die "unable to open grep pipe:";
    while (my $line = <GR>) {
      chomp $line;
      $line =~ /(.*)/;
      print STATS "num_pfam_hits\t$1\n";
      }
  } else {
    print "ERROR: either $data_file does not exist or is empty!";
    print STDERR "ERROR: either $data_file does not exist or is empty!";
    exit 0;
  }

  $cmd = 'grep PF ' . $data_file . ' | cut -f1,1 -d " " | sort -u | wc -l';
  open (GR, "$cmd |") or die "unable to open grep pipe:";
  while (my $line = <GR>) {
    chomp $line;
    $line =~ /(.*)/;
    print STATS "num_transcripts_with_pfam_hits\t$1\n";
  }
} elsif ($job_type eq $QUANTIFICATION) {
  $data_file = $project_name . $QUANTIFICATION_EXT;
  if ((-e $data_file) && (not -z $data_file)) {
    print STATS "quantification_done\t'Y'\n";    
  } else {
    print STATS "quantification_done\t'N'\n";
  }
} else {
  print "ERROR: UNKNOW JOB_TYPE: $job_type\n";
  print STDERR "ERROR: UNKNOW JOB_TYPE: $job_type\n";
  exit 0;

}

close GR;
close STATS
#
#c10node12.acis.ufl.edu{pwilliams}637: cat *.stats num_assembled_seqs
#7592 num_assembled_bases 3232991 num_transcripts_with_blast_nr_hits
#-547858 num_blast_nr_hits 20977 c10node12.acis.ufl.edu{pwilliams}638:
#grep Query= *nr.txt | wc -l 7440 c10node12.acis.ufl.edu{pwilliams}639:
#grep '>' *nr.txt | sort -u | wc -l 15586
#c10node12.acis.ufl.edu{pwilliams}640: grep '\*\*\*\*' *nr.txt | wc -l
#2095 =cut

#  55018:37cut -f1,1 *panther* | sort -u | wc -l
# 55118:37cut -f1,1 *panther* | wc -l
# 55218:38wc -l *pant*
