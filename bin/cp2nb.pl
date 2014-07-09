#!/usr/bin/perl -w
use strict;
use File::Copy;

my ($proj,$nb) = @ARGV;
if ((not defined $proj) || (not defined $nb)) {
    print "\nUsage: $0 <proj> <nb_num: 2 (hpc), 3(oem), 4(pleurobrachia), or 5(aplysia)>\n";
    exit 0;
}

if (($nb ne 2)&& ($nb ne 3)&& ($nb ne 4)&& ($nb ne 5)) {
   print "invalid <nb_num>: $nb\n";
   exit 0;
}

if ($nb == 2) { $nb = 'hpc2'; }
elsif ($nb == 3) { $nb = 'oem'; }
elsif ($nb == 4) { $nb = 'pb'; }
elsif ($nb == 5) { $nb = 'ap'; }
else { print "$nb is not in oem, hpc2, ap, pb\n"; exit 0; }

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

$cmd = "/home/pwilliams/bin/rsy.pub $proj $nb";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "rm -fr $proj";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
