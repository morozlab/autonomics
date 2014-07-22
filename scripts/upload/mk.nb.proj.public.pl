#!/usr/bin/perl -w
use strict;

my ($pid) = @ARGV;
if (not defined $pid) {
print "\nUsage: $0 pid\n";
exit 0;
}

my $tend =  " | mysql -u moroz_lab --password=Whitney2011 --database=moroz_lab";

my $cmd = "echo \"select project_name from project_directory where projectID = $pid\"" . $tend;
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "echo \"delete from permissions where project_id = $pid\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

$cmd = "echo \"insert into permissions values (2, $pid)\"" . $tend;
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
