#!/usr/bin/perl
use strict;
use warnings;

my $cmd;

my @projs = qw(

hlc4

);

foreach my $p (@projs) {
 $cmd = "loadnb.pl $p NT ALL";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}
