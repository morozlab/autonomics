#!/usr/bin/perl
use strict;
use warnings;

my $MAX_NB = 5;

my ($p,$ntaa,$nb) = @ARGV;
if ((not defined $p) || (not defined $ntaa) || (not defined $nb)) {
   print "\nUsage: $0 proj_name <NT or AA> <NB#>\n";
   exit 0;
}

if ((not ($ntaa == 'NT')) && (not ($ntaa == 'AA'))) {
   print "\nUsage: $0 proj_name <NT or AA> <NB#>\n";
   exit 0;
}

$nb = int $nb;
if (($nb < 1) || ($nb > $MAX_NB)) {
   print "NB# must be an int between 1 and $MAX_NB\";
   exit 0;
}

if (($ntaa == 'NT')) && (not ($ntaa == 'AA'))) {
   print "\nUsage: $0 proj_name <NT or AA> <NB#>\n";
   exit 0;
}

my $cmd = "load_project.pl $p $ntaa $nb";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
