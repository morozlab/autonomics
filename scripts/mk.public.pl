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
41
5
87
11
84
9
83
4
88
94
86
8
89
27
39
92
95
96
12
97
98
14
99
100
74
78
79
75
76
80
15
81
40
91
26
38
6
85
);

# foreach my $pid (41,5,87,11,84,9,83,4,88,3,86,8,89,27,39,92,19,13,12,17,20,14,18,21,74,78,79,80,38,15,75,76,26,81,40,91,85,6) {

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
