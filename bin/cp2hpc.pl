#!/usr/bin/perl -w
use strict;
use File::Copy;

my ($proj,$un) = @ARGV;
if ((not defined $proj)|| (not defined $un)) {
    print "\nUsage: $0 <proj> <your_gator_user_name>\n";
    exit 0;
}

if($proj =~ m/.*\/$/) {
    $proj =~ s/\/$//;
}

$proj =~ s/\///;

my $dir = "/srv/data2/pipeline/";
my $dir2 = $dir . $proj;
chdir $dir2 or die "unable to chdir $dir2";

if (not -e $proj) {
  unless(mkdir $proj) { die "Unable to create directory: $proj"; }
}

my $p = $proj . "_project.fasta";
my $cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_pfam.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_blast_nr.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_blast_swissprot.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_KEGG.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_GO.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_gocats.txt";
$cmd = "rsync $p $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$p = $proj . "_quantification.txt";
if (-e $p ) {
  $cmd = "rsync $p $proj";
  print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

$cmd = "rsync -avzl $proj $un\@hipergator.hpc.ufl.edu:/rlts/moroz/data2load";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "rm -fr $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
