#!/usr/bin/perl
use strict;
use warnings;

my $apath  = $ENV{AUTONOMICS_PATH };
my $config_file = $apath . '/config.pl';
my %config = do $config_file;
my $pw = $config{db_root_password};

my $str = 'delete from load_info where project_id = <NN> and data = quant;';


my ($proj_name, $data_type, $nb) = @ARGV;
if( (not defined $proj_name) ||
    (not defined $nb) ||
    (not defined $data_type)) {
    print "\nUsage: $0 <proj_name> <NT || AA> <NB_num>\n";
    exit 0;
}

print "\n\nif (bad) quantif data already loaded, you need to first do in mysql:\n\n\t$str\n\n";

exit 0;

if (($data_type ne "AA") && ($data_type ne "NT")) {
  print "data_type must be AA or NT, you entered: $data_type\n";
  exit 0;
}

my $cmd = "loadnb.pl $proj_name $data_type QUANT $pw $nb";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
