#!/usr/bin/perl
use strict;
use warnings;

my $cmd;

my @projs = qw(

Capitella_proteins
Helobdella_proteins

);

foreach my $p (@projs) {
 $cmd = "load.project.pl $p AA ALL";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}
