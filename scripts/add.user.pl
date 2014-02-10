#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my ($name, $pw, $email) = @ARGV;
if ((not defined $name) ||
    (not defined $pw) ||
(not defined $email)     ) {
    print "\nUsage: $0 uname  pw email\n";
    exit 0;
}

my $tend =  " | mysql -u moroz_lab --password=Whitney2011 --database=moroz_lab";

 my $cmd = "echo \"insert into users values(DEFAULT, \'$name\', \'$pw\', \'$name\', \'$email\')\"" . $tend;
 print "$cmd\n";
 system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }
# hawkins 7272ca19f8f236951c382d87ae2cfac1 
