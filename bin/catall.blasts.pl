#!/usr/bin/perl -w
use strict;



my ($projname, $nrsw, $num_last) = @ARGV;
if ((not defined $projname) ||
    (not defined $nrsw) ||
    (not defined $num_last)) {

    print "\nUsage: $0 <proj_name> <NR | SW | PFAM> <num_last_out_file>\n";
    exit 0;
}

my $pfile = "";
my $ending = "";
my $outfile = "";
if ($nrsw eq 'NR') {
  $pfile = "_blast_nr";
  $ending = $pfile . ".txt";
  $outfile =  $projname . $ending;
} elsif  ($nrsw eq 'SW') {
  $pfile = "_blast_swissprot";
  $ending = $pfile . ".txt";
  $outfile =  $projname . $ending;
} elsif  ($nrsw eq 'PFAM') {
  $pfile = "_pfam";
  $ending = ".stdout";
  $outfile =  $projname . "_pfam.txt";
} else {
    print "2nd arg must be NR or SW or PFAM\n";
  exit 0
}


if (-e $outfile) {
    print "$outfile already exists\n";
    exit 0;
}

my $cmd = "touch $outfile";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

foreach my $num (1 .. $num_last) {
  my $infile =  $projname . $pfile . $num . $ending;
  my $cmd = "cat $infile >> $outfile";
  print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

