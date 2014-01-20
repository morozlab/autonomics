#!/usr/bin/perl
use strict;
use warnings;

my $cmd;

my @projs = qw(

Aplysia_LE_single_cell_005_clcbio
);

foreach my $p (@projs) {
 $cmd = "load.project.pl $p NT ALL";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}
