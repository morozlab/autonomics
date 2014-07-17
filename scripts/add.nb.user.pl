#!/usr/bin/env perl
use strict;
use nb_names;

my ($uname, $rname, $pw, $email, $nb_abrev) = @ARGV;
if ((not defined $uname) ||
    (not defined $pw) ||
    (not defined $nb_abrev) ||
    (not defined $rname) ||
#    (not defined $pid) ||
    (not defined $email)) {
  print "\nUsage: $0  <uname>  <rname>  <pw>  <email> <nb_abrev> \n";
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

my $tend =  " | mysql -u nb --password=q8yqJ6zk -h db1.ufhpc --database=$nb_name";

my $cmd = "echo \"insert into users values(DEFAULT, \'$uname\', \'$pw\',\'$rname\', \'$email\')\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "perl /rlts/moroz/autonomics/bin/sqn  $nb_abrev  uid  $uname  | tail -n 1 | cut -f1,1";
open (GR, "$cmd |") or die "unable to open  pipe:";
my $uid = <GR>; chomp $uid;
close GR;

