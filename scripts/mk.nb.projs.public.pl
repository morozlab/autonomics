#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;

my ($file) = @ARGV;
if ((not defined $file)) {
    print "\nUsage: $0 <filename-with-list-of-projects to load into nb\n";
    exit 0;
}

if (not -e $file) {
    print "file: $file does not exist\n";
  exit 0;
}

my $cmd = "foreach.line.in.file.do.cmd $file mk.nb.proj.public.pl";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }


