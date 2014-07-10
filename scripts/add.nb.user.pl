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


my ($uname, $rname, $pw, $email, $pid, $nb_abrev) = @ARGV;
if ((not defined $uname) ||
    (not defined $pw) ||
    (not defined $nb_abrev) ||
    (not defined $rname) ||
    (not defined $pid) ||
    (not defined $email)) {
  print "\nUsage: $0  <uname>  <rname>  <pw>  <email>  <PID or -1 if GLOBAL>  <nb_abrev> \n";
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


my $cmd = "echo \"insert into users values(DEFAULT, \'$uname\', \'$pw\',\'$rname\', \'$email\')\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "perl /rlts/moroz/autonomics/bin/sqn  $nb_abrev  uid  $uname  | tail -n 1 | cut -f1,1";
open (GR, "$cmd |") or die "unable to open  pipe:";
my $uid = <GR>; chomp $uid;
close GR;

if ($pid == -1) {
  $cmd = "echo \"delete from permissions where user_id = \'$uid\'\"" . $tend;
  print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

$cmd = "echo \"insert into permissions values(\'$uid\', \'$pid\')\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
