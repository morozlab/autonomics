#!/usr/bin/perl
use strict;
use warnings;

my ($pname, $seq_type, $load_type, $debug) = @ARGV;
if ((not defined $pname) ||
    (not defined $seq_type) ||
    (not defined $load_type)) {
    print "\nUsage: $0  <proj_name>  <seq_type[NT|AA]>  <load_type[BASIC|NR|SP|GO|GOCAT|KEGG|BESTA]>   <DEBUG>\n";
    exit 0;
}

if (not defined $debug) { $debug = ""; }
else { 
  if ($debug ne 'DEBUG') {
    print "<DEBUG> must be DEBUG, you used: $debug\n";
  exit 0;
  }
  $debug = "--debug";
}

$seq_type = uc($seq_type);
$load_type = uc($load_type);

if (($seq_type ne 'AA') && ($seq_type ne 'NT')) {
  print "<seq_type> must be AA or NT, you used: $seq_type\n";
  exit 0;
}
if (($load_type ne 'BASIC') && ($load_type ne 'NR') && ($load_type ne 'SP') && ($load_type ne 'GO')  && ($load_type ne 'KEGG')   && ($load_type ne 'BESTA')  && ($load_type ne 'GOCAT')) {
  print "<load_type> must be BASIC or NR or SP or GO or GOCAT or KEGG or BESTA, you used: $load_type\n";
  exit 0;
}

my $host  = $ENV{HOST};

my @date = `date`;
chomp(@date);
my $pwd = `pwd`;
chomp($pwd);
print "Executed as \"$0 @ARGV\" on $ENV{HOST} in $pwd  @date\n";

if ($host =~ /hpc/) { $host = 'hpc'; }
elsif ($host =~ /oem/) { $host = 'oem'; }
else {  print "unable to find host from $host\n"; exit 0; }

my $dir = "";
my $sdir = "";

if ($host eq 'oem') {
  $dir = '/data/neurobase_load_data/';
  $sdir = '/home/oem/python_software/autonomics/scripts';
}
if ($host eq 'hpc') {
  $dir = '/data/autonomics_data/';
  $sdir = '/home/pwilliams/python_software/autonomics/scripts';
}

if ($pwd ne $sdir) { print "You need to be in $sdir to run this script\n"; exit 0; }

my $log = "log." . $pname . "." . $load_type;
if (-e $log) {  print "log file: $log already exists\n"; exit 0; }

if ($seq_type eq 'AA') { $seq_type = '--seq-type AA'; }
else { $seq_type = ""; }

my $quant;
my $quant_file = $dir . "/" . $pname . "/" . $pname . "_quantification.txt";
print "q: $quant_file";

if (-e $quant_file) {
    $quant = "";
}  else {
    $quant = "--no-quant";
}

if ($load_type eq 'NR') { $load_type = '--nr --annotation'; }
elsif ($load_type eq 'SP') { $load_type = '--swissprot --annotation'; }
elsif ($load_type eq 'GO') { $load_type = '--GO'; }
elsif ($load_type eq 'KEGG') { $load_type = '--KEGG'; }
elsif ($load_type eq 'BEST') { $load_type = '--best-annotation'; }
elsif ($load_type eq 'GOCAT') { $load_type = '--GOCategories'; }
else { $load_type = ""; }

my $cmd = "python upload.py -p $pname --user root --password Moof2011 -v --data-dir $dir $seq_type $load_type $debug $quant >& $log &";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }
