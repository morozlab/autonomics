#!/usr/bin/perl -w
use strict;

my ($f,$start,$end,$nrsp) = @ARGV;
if ((not defined $f) || (not defined $start)  || (not defined $end) || (not defined $nrsp)) {
    print "\nUsage: $0 <blast_output_file (leave off NNN & what follows> <start_num> <end_num> <nr | sp>\n";
    exit 0;
}

my $suffix;

if ($nrsp eq 'nr') {
  $suffix = "_blast_nr.txt";
} elsif ($nrsp eq 'sp') {
 $suffix = "_blast_swissprot.txt";
} else {
  print  "last arg must be nr or sp\n";
  exit 0;
}

foreach my $num ($start .. $end) {
    my $file = $f . $num . $suffix;
    my $res = `grep Query= $file | wc -l`;
    chomp $res;
    print "$file\t$res\n";
}

