#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my ($x) = @ARGV;
if ((not defined $x)) {
    print "\nUsage: $0 1\n";
    exit 0;
}

my $tend =  " | mysql -u moroz_lab --password=Whitney2011 --database=moroz_lab";

foreach my $pid (41,5,11,25,9,4,3,8,27,39,19,13,12,17,20,14,18,21,74,78,79,80,38,15,75,76,26,81,40,6,83,84,85,86) {

 my $cmd = "echo \"select project_name from project_directory where projectID = $pid\"" . $tend;
## print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

 $cmd = "echo \"delete from permissions where project_id = $pid\"" . $tend;
## print "$cmd\n";
 #system($cmd);
 #if ( $? ) { die "Command failed: $cmd: $!"; }

 $cmd = "echo \"insert into permissions values (2, $pid)\"" . $tend;
## print "$cmd\n";
 #system($cmd);
 #if ( $? ) { die "Command failed: $cmd: $!"; }
}
