#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my ($user_id, $pid) = @ARGV;
if ((not defined $user_id) ||
    (not defined $pid) ) {
    print "\nUsage: $0 user_id  pid\n";
    exit 0;
}

my $tend =  " | mysql -u moroz_lab --password=Whitney2011 --database=moroz_lab";

 my $cmd = "echo \"insert into permissions values(\'$user_id\', \'$pid\')\"" . $tend;
 print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }
