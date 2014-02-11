#!/usr/bin/perl
use strict;
use warnings;

my $cmd;

my ($p,$ntaa) = @ARGV;
if ((not defined $p) || (not defined $ntaa)) {
   print "\nUsage: $0 proj_name <NT or AA>\n";
   exit 0;
}

$cmd = "load_project.pl $p $ntaa ALL";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
