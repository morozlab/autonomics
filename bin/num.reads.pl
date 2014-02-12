#!/usr/bin/perl -w
use strict;

my ($fq) = @ARGV;
if ((not defined $fq)) {
    print "\nUsage: $0 <fastq file>\n";
    print " You need to be in dir where these files are \n";
    exit 0;
}

my $cmd = "wc -l $fq";
my $res =`$cmd`;
my( $r, $x ) = split( /\s+/, $res );
$r = $r /4;
print "$r\n";
