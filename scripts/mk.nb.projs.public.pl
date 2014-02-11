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

my @pids = qw(


);

foreach my $pid (@pids) {

 my $cmd = "echo \"select project_name from project_directory where projectID = $pid\"" . $tend;
## print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

 $cmd = "echo \"delete from permissions where project_id = $pid\"" . $tend;
## print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }

 $cmd = "echo \"insert into permissions values (2, $pid)\"" . $tend;
## print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }
}
