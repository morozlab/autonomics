#!/usr/bin/perl
use strict;
use warnings;

my ($pname, $seq_type, $load_type, $rootpw, $nb, $debug) = @ARGV;
if ((not defined $pname) ||
    (not defined $seq_type) ||
    (not defined $nb) ||
    (not defined $rootpw) ||
    (not defined $load_type)) {
    print "\nUsage: $0  <proj_name>  <seq_type[NT|AA]>  <load_type[ALL|BASIC|NR|ALIGN_NR|ALIGN_SP|SP|PFAM|GO|GOCAT|KEGG|BESTA|QUANT]> <root_password> <NB#>  <DEBUG --optional>\n";
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
if (($load_type ne 'BASIC') && ($load_type ne 'NR')&& ($load_type ne 'ALL') && ($load_type ne 'SP') && ($load_type ne 'GO') && ($load_type ne 'QUANT') && ($load_type ne 'ALIGN_NR')&& ($load_type ne 'ALIGN_SP')  && ($load_type ne 'KEGG')   && ($load_type ne 'BESTA')   && ($load_type ne 'PFAM')  && ($load_type ne 'GOCAT')) {
  print "<load_type> must be BASIC or NR or ALL or ALIGN_NR or ALIGN_SP or SP or GO or GOCAT or KEGG or BESTA or PFAM, you used: $load_type\n";
  exit 0;
}

my @date = `date`;
chomp(@date);
my $pwd = `pwd`;
chomp($pwd);

my $dir = $ENV{NEUROBASE_LOAD_DATA_PATH};
my $sdir = $ENV{AUTONOMICS_PATH} . '/scripts';

if (not -e "upload.py") { print "You need to be in $sdir to run this script\n"; exit 0; }

my $log = "log." . $pname . "." . $load_type;
if (-e $log) {  print "log file: $log already exists\n"; exit 0; }

if ($seq_type eq 'AA') { $seq_type = '--seq-type AA'; }
else { $seq_type = ""; }

my $quant;
my $quant_file = $dir . "/" . $pname . "/" . $pname . "_quantification.txt";

if (-e $quant_file) {
    $quant = "";
}  else {
    $quant = "--no-quant";
}

print "load_type: $load_type\n";

if ($load_type ne 'ALL') {
  if ($load_type eq 'NR') { $load_type = '--nr --annotation'; }
  elsif ($load_type eq 'ALIGN_NR') { $load_type = '--nr --alignments'; }
  elsif ($load_type eq 'ALIGN_SP') { $load_type = '--swissprot --alignments'; }
  elsif ($load_type eq 'SP') { $load_type = '--swissprot --annotation'; }
  elsif ($load_type eq 'GO') { $load_type = '--GO'; }
  elsif ($load_type eq 'PFAM') { $load_type = '--pfam'; }
  elsif ($load_type eq 'QUANT') { $load_type = '--quantification'; }
  elsif ($load_type eq 'KEGG') { $load_type = '--KEGG'; }
  elsif ($load_type eq 'BEST') { $load_type = '--best-annotation'; }
  elsif ($load_type eq 'GOCAT') { $load_type = '--GOCategories'; }
  else { $load_type = ""; }

  my $cmd = "python upload.py -p $pname --user root --password $rootpw -v $seq_type $load_type --nb $nb  $debug $quant >& $log";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

} else {

  my $log = "log." . $pname . ".BASIC";
  $load_type = "";
  my $cmd = "python upload.py -p $pname --user root --password $rootpw -v  $seq_type $load_type --nb $nb  $debug $quant >& $log";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

  print "done #1 starting #2\n";

  $load_type = '--nr --annotation';
  $log = "log." . $pname . ".NR";
  $cmd = "python upload.py -p $pname --user root --password $rootpw -v $seq_type $load_type --nb $nb  $debug $quant >& $log";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

  print "done #2 starting #3\n";

  $log = "log." . $pname . ".ALIGN_NR";
  $load_type = '--nr --alignments';
  $cmd = "python upload.py -p $pname --user root --password $rootpw -v $seq_type $load_type --nb $nb  $debug $quant >& $log";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

  $cmd = "add2gp $pname";
  print "$cmd\n";
  system ($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }

}
