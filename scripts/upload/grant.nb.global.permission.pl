#!/usr/bin/env perl
use strict;
use nb_names;
my ($user_id,$nb_abrev) = @ARGV;
if ((not defined $user_id) ||
    (not defined $nb_abrev)) {
  print "\nUsage: $0 user_id nb_abrev \n";
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

my $cmd = "echo \"delete from permissions where user_id = \'$user_id\'\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "echo \"insert into permissions values(\'$user_id\', -1)\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

