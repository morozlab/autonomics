#!/usr/bin/env perl
use strict;
use warnings;

my %nb_abrevs = (
'misc' => 'misc',
'apt' => 'aplysia',
'apn' => 'aplysia2',
'pb' => 'pleurobrachia',
'sandbox' => 'sandbox',
'genomes' => 'genomes',
'sponges' => 'porifera',
'nb1' => 'nb1',
'molluscs' => 'molluscs',
#'snails' => 'gastropoda',
'secr' => 'secretoryMolecules',
'verts' => 'vertebrates',
'squid' => 'cephalopods',
'combs' => 'ctenophora',
    );

my ($p,$ntaa,$nb_abrev) = @ARGV;
if ((not defined $p) || (not defined $ntaa) || (not defined $nb_abrev)) {
   print "\nUsage: $0 proj_name <NT or AA> <NB_abrev, e.g. pb>\n";
   print "NB_abrev  NB_name\n";
   while( my( $key, $value ) = each %nb_abrevs ){
	print "$key\t $value\n";
   }
   exit 0;
}

my $nb_name;
if( exists($nb_abrevs{$nb_abrev} ) ){
    $nb_name = "nb_".$nb_abrevs{$nb_abrev};
    print "nb_name: $nb_name\n";
} else {
    print "\nUsage: $0 proj_name <NT or AA> <NB_abrev, e.g. pb>\n";
    print "\nINVALID NB_ABREV: $nb_abrev    here are the valid ones\n\n";
    print "NB_abrev  NB_name\n";
    while( my( $key, $value ) = each %nb_abrevs ){
	print "$key\t $value\n";
    }
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
