#!/usr/bin/env perl
use strict;
use warnings;
use nb_names;

# /rlts/moroz/autonomics/scripts/

my ($p,$ntaa,$nb_abrev) = @ARGV;
if ((not defined $p) || (not defined $ntaa) || (not defined $nb_abrev)) {
   print "\nUsage: $0 proj_name <NT or AA> <NB_abrev, e.g. pb>\n";
   print "NB_abrev  NB_name\n";
   nb_names::print_valid_nb_abrevs();
   exit 0;
}

my $nb_name = nb_names::nb_abrev_lookup($nb_abrev);

if ($nb_name) {
    print "nb_name: $nb_name\n";
} else {
    print "\nUsage: $0 proj_name <NT or AA> <NB_abrev, e.g. pb>\n";
    print "\nINVALID NB_ABREV: $nb_abrev    here are the valid ones\n\n";
    print "NB_abrev  NB_name\n";
    nb_names::print_valid_nb_abrevs();
    print "\n";
    exit 0;
}

if (($ntaa ne 'NT') && ($ntaa ne 'AA')) {
   print "\nUsage: $0 proj_name <NT or AA> <NB_num>\n";
   exit 0;
}

my $cmd = "load_project.pl $p $ntaa $nb_name $nb_abrev";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
