#!/usr/bin/env perl
use strict;

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


my ($nb_abrev, $pid) = @ARGV;
if ((not defined $pid) ||
    (not defined $nb_abrev)) {
  print "\nUsage: $0 nb_abrev pid\n";
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

my $tend =  " | mysql -u nb --password=q8yqJ6zk -h db1.ufhpc --database=$nb_name";

my $cmd = "echo \"delete from permissions where project_id = \'$pid\'\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "echo \"insert into permissions values(2, \'$pid\')\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

