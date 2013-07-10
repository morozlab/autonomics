#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my($proj, $mira, $data, $paired_end, $noass);

my $cmd_line = "$0 @ARGV";
my $args0 = "@ARGV";

GetOptions(
  "proj=s" => \$proj,
  "mira" => \$mira,
  "data=s" => \$data,
  "noass" => \$noass,
  "paired_end" => \$paired_end,
);

sub printUsage{
  print "$0\n";
  print "usage:\n";
  print "	-proj :project_name (REQUIRED)\n";
  print "       -mira   (OPTIONAL)\n";
  print "       -paired_end   (OPTIONAL)\n";
  print "       -noass   (OPTIONAL)\n";
  print "       -data   (OPTIONAL [ NT | AA ] - REQUIRED if -noass used)\n";
};

if(!defined($proj)) { printUsage(); exit 0; }
if(!defined($mira)) { $mira = 0; }  else { $mira = 1; }
if(!defined($noass)) { $noass = 0; }  else { $noass = 1; }
if(!defined($paired_end)) { $paired_end = 0; } else { $paired_end = 1; }

if ($noass && (!defined($data))) {
  print "-data NT | AA must be specified when -noass is set\n";
  exit 0;
}

if ($noass && (($data ne "NT") && ($data ne "AA"))){
  print "-data must be NT or AA\n";
  exit 0;
}

my $blastp;

if ($noass) {
  if ($data eq "AA") { $blastp = 1; }
  else { $blastp = 0; }
}

if ($paired_end && $mira) {
  print " can not set paired_end and mira flags at same time\n";
  exit 0;
}

if (($paired_end || $mira) && $noass) {
  print " can not set paired_end or mira flags with noass flag\n";
  exit 0;
}

my $arg = "";
if ($mira) {
  $arg = "--set-args \"assemble|pipeline_args;--assembler mira\"";
}

my $paired_end_arg = "";
if ($paired_end) {
  $arg = "--set-args \"quality_trim|pipeline_args;--paired-end\" \"adapter_trim|pipeline_args;--paired-end\" \"read_normalization|pipeline_args;--paired-end\" \"assemble|pipeline_args;--paired-end\"";
}

my $fa;
my $fa2;
my $projdir = "/srv/data2/pipeline/" . $proj . "/";

if (!-d $projdir) {
  print "$projdir does not exist\n";
  exit 0;
}

if ($noass) {
  $fa = $projdir . $proj . "_project.fasta";
} else {
  $fa = $projdir . $proj . "/" . $proj . ".fastq";
  $fa2 = $projdir . $proj . "/" . $proj . ".fastq.end2";
  if ((-e $fa2) && (!$paired_end)) {
    print "this project has paired end data but the paired_end flag is not set!\n";
    exit 0;
  }
}

if (not -e $fa) { 
  print "$fa is REQUIRED but does not exist!\n";
  exit 0;
}

if ($paired_end) {
  if (not -e $fa2) {
    print "$fa2 is REQUIRED but does not exist!\n";
    exit 0;
  }
}

if ($blastp) {
  $arg = "--set-args \"blast_nr|pipeline_args;--aligner blastp\" \"blast_swissprot|pipeline_args;--aligner blastp\" \"pfam|pipeline_args;\"";
}

my $cmd;
if ($noass) {
  $cmd = "python /home/pwilliams/python_software/autonomics/scripts/systemtools.py --add-project --assign-workflow --add-jobs blast_nr blast_swissprot pfam kegg go --set-config blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+ $arg --project-names $proj";
} else {
  $cmd = "python /home/pwilliams/python_software/autonomics/scripts/systemtools.py --add-project --assign-workflow --add-jobs adapter_trim quality_trim read_normalization assemble quantification blast_nr blast_swissprot pfam kegg go --set-config adapter_trim:+ quality_trim:+ read_normalization:+ assemble:+ quantification:+ blast_nr:+ blast_swissprot:+ pfam:+ kegg:+ go:+ $arg --project-names $proj";
}

print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

print "job submitted successfully\n";

