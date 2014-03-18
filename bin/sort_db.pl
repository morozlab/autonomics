#!/usr/bin/perl -w

use strict;
use IO::File;
use File::Basename;

#####################################################################
#
#                    sort_db.pl
#                    ==========
#
#           Peter Williams    LLNL    May 2007
#
# Using an out-of-core bin sort, sequences in a very large fasta
# database are sorted in ascending order of sequence length.  To sort
# a 20GB database, takes about 45 minutes on an AMD node.
#
# sort_db.pl <database_with_path> <output_dir_with_path>
#
# The sorted database will be written to:
#     <output_dir_with_path>/<database_name>/<database_name>_sorted

# For example:
#   sort_db.pl /usr/db/NT /usr/local/partitioned_databases
#      will create the dir /usr/local/partitioned_databases/NT and
#      then write out /usr/local/partitioned_databases/NT/NT_sorted
#
# This particular directory structure is important to subsequent steps
# such as partition_db.pl and the later use of the partitioned databases.
# For more information, see the file README_pBLAST.
#
# NB: the output directory MUST be a different directory than the one
#     in which the database lives. 
#     So:  sort_db.pl /usr/db/NT /usr/db is NOT allowed.
#
# Methodology:
#   chooses bin limits (see next paragraph) based on distribution
#   of sequence lengths, such that each bin is less than $MAX_BIN_SIZE
#
#   sequences are then placed in bins (a bin is a disk file).  So,
#   e.g. if the bin limits are (5_000 30_000 100_000) then there are 4
#   bins:
#     bin0 has all sequences with lengths from 0 to 5_000
#     bin1 has all sequences with lengths from 5_001 to 30_000
#     bin2 has sequences with lengths from 30_001 to 100_000
#     bin3 has sequences with lengths from 100_001 and up.
#
#  then each bin, starting with bin0, is sorted in-core, and results
#  concatenated to the database_sorted disk file.
#
#####################################################################


#####################################################################
#
#   Set max bin size --- 2GB is good if you have at least 4GB of RAM.
#   otherwise set it to about 50% of the amount of your RAM.
#
# Note: MAX_BIN_SIZE must be greater than the length of longest 
# sequence in the database.  You can use "scandb <database>" to find 
# the length of longest sequence.
#
#####################################################################

my $MAX_BIN_SIZE = 2 * 1_073_741_824; # 2GB



my ($db_file_with_path, $output_dir) = @ARGV;
if ((not defined $db_file_with_path) ||
    (not defined $output_dir)) {
  print "\nUsage: $0 <database_with_path> <output_dir_with_path>\n\n";
  print "e.g. $0 /usr/db/NT /usr/local/partitioned_databases\n";
  print "\nwill create the dir /usr/local/partitioned_databases/NT and\n";
  print "then write out /usr/local/partitioned_databases/NT/NT_sorted\n\n";
  exit;
}

# print "MAX_BIN_SIZE = $MAX_BIN_SIZE\n";

my $db_name = basename($db_file_with_path);
my $db_dir =  dirname($db_file_with_path);
if ($output_dir =~ /(.*)\/$/) { # remove trailing "/" if any.
  $output_dir = $1
}

if ($db_dir eq $output_dir) {
#    die "\nError:output_directory must be different than the directory in which $db_name lives\n";
}

# $output_dir .= '/' . $db_name;
$output_dir = $db_dir;

# unless (-e $output_dir) { print `mkdir $output_dir`; }

my $temp_dir = $output_dir . '/temp';
my $sorted_db_file = $output_dir . '/' . $db_name . '_sorted';

unless (-e $temp_dir) { print `mkdir $temp_dir`; }

#####################################################################
#
#     get length of each sequence => @sizes
#
#####################################################################

my $sequence = undef;
my $id;
my $next_id = "";
my $line;
my @sizes;

# print "Starting Phase 1 ...\n";
my $start_time = time;

if (-d $db_file_with_path) { die "\nCan't open $db_file_with_path -- it is a directory\n\n"; }
else { open INDATA,"$db_file_with_path" or die "Can't open $db_file_with_path: $!\n"; }

while (1) {
  $line = <INDATA>;
  if ((defined $line) && ($line =~ /^\s+$/)) { next; } # skip blank lines
  if ((!defined($line)) || ($line =~ /^(>.*)/)) {
    $id = $next_id;
    if (defined $line) { $next_id = $1; }
    if (defined $sequence) { # sequence is undef at start
      push @sizes, (length $sequence);
    }
    $sequence = "";
    if (!defined($line)) {
      last;
    }
  } else {
    $sequence .= $line;
  }
}

close INDATA;

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
$start_time = time;

######################################################################
#
# sort @sizes => @sorted_sizes
#
######################################################################
# print "Starting Phase 2 ...\n";

my @sorted_sizes = sort {$a <=> $b} @sizes;

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
$start_time = time;

######################################################################
#
# Find bin limits by scanning @sorted_sizes so each bin, when filled,
# will be < $MAX_BIN_SIZE
#
# sequences size:
# 0           limits[0], limits[0]+1 limits[1],  .. limits[$#limits]  infinity
# |--db_$limits[0]----|  |---db_$limits[1]---|   .. |----------db_inf--------|
#
# if @limits = (5_000 30_000) then there are 3 bins:
# bin0 has sequences with lengths from 0 to 5_000
# bin1 has sequences with lengths from 5_001 to 30_000
# bin2 has sequences with lengths from 30_001 and up
#
# Assumption: no set of equal length seqs has cardinality > MAX_BIN_SIZE.
#
######################################################################
# print "Starting Phase 3 ...\n";

my $starting_a_bin = 1;
my $last_size = 0;
my $last_diff_size;
my $sum = 0;
my $count = 0;
my @limits;

while (my $seq_size = shift @sorted_sizes) {

  # next logic required to avoid overflow when encounter multiple seqs
  # of equal length.
  if ($seq_size != $last_size) { $last_diff_size = $last_size; $count=1; }
  else { $count++; }

  if ($starting_a_bin) { $starting_a_bin = 0; }
  $sum += $seq_size;
  if ($sum > $MAX_BIN_SIZE) {
    $starting_a_bin = 1;
    if ($seq_size != $last_size) {
      push @limits, $last_size;
      $sum = $seq_size;
    } else { # encountered several seqs of equal len
      push @limits, $last_diff_size;
      $sum = $seq_size * $count;
    }
  }
  $last_size = $seq_size;
}

my $num_limits = scalar @limits;

 if ($num_limits > 0) { print "limits: @limits\n"; }
 else { print "limits: infinity\n"; }

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
$start_time = time;

######################################################################
#
# open file for each bin
#
######################################################################
# print "Starting Phase 4 ...\n";

my @fha;
for (my $idx=0; $idx<$num_limits; $idx++) {
  my $file_name = $temp_dir . '/' . $db_name . '_' . $limits[$idx];
  push @fha, new IO::File ">$file_name";
}

my $file_name = $temp_dir . '/' . $db_name . '_inf';
push @fha, new IO::File ">$file_name";

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
$start_time = time;

######################################################################
#
# distribute sequences in database to the bins
#
######################################################################

# print "Starting Phase 5 ...\n";

open INDATA,"$db_file_with_path" or die "Can't open $db_file_with_path: $!\n";

$sequence = undef;
$next_id = "";

my $genome_num = 0;

while (1) {
  $line = <INDATA>;
  if ((defined $line) && ($line =~ /^\s+$/)) { next; }
  if ((!defined($line)) || ($line =~ /^(>.*)/)) {
   $id = $next_id;
   if (defined $line) { $next_id = $1; chomp $next_id; }
    if (defined $sequence) {
     &distribute_seq($id, uc($sequence), $genome_num);
     $genome_num++;
   }
   $sequence = "";
   if (!defined($line)) {
     last;
   }
  } else {
    $sequence .= $line;
  }
}

close INDATA;
for (my $idx=0; $idx<=$#fha; $idx++) {
  close $fha[$idx];
}

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
$start_time = time;

######################################################################
#
# sort seqs in each bin & concatenate sorted seqs to $db_file_with_path."_sorted"
#
######################################################################
# print "Starting Phase 6 ...\n";

my @seqs;
my %hash;

open OUTDATA,">$sorted_db_file" or die "Can't open $sorted_db_file: $!\n";

@fha=();

for (my $idx=0; $idx<$num_limits; $idx++) {
  my $file_name = $temp_dir . '/' . $db_name . '_' . $limits[$idx];
  push @fha, new IO::File "$file_name";
}
$file_name = $temp_dir . '/' . $db_name . '_inf';
push @fha, new IO::File "$file_name";

for (my $idx=0; $idx<=$#fha; $idx++) {
  my $sequence;
  my $id;
  my $next_id = "";
  my $genome_num=0;
  my $line;

  @seqs = ();
  %hash = ();

  while (1) {
    $line = readline($fha[$idx]);
    if ((defined $line) && ($line =~ /^\s+$/)) { next; }
    if ((!defined($line)) || ($line =~ /^(>.*)/)) {
      $id = $next_id;
      if (defined $line) { $next_id = $1; chomp $next_id; }
      if (defined $sequence) {
        $seqs[$genome_num] =  "$id\n$sequence";
        $hash{$genome_num} = length $sequence;
	$genome_num++;
      }
      $sequence = "";
      if (!defined($line)) {
	last;
      }
    } else {
      $sequence .= $line;
    }
  }

  foreach my $key (sort {$hash{$a} <=> $hash{$b} } keys %hash) {
    print OUTDATA "$seqs[$key]";
  }
  close $fha[$idx];
}

# printf  "                 Completed in: %8.3f mins.\n", ((time - $start_time) / 60.0);
 my $cmdx = "/bin/rm -fr $temp_dir";
 print "Removing temp dir..........";
 system($cmdx);
 if ( $? ) { die "Command failed: $cmdx: $!"; }
 print "done\n";


######################################################################
#
#                    subroutines
#
######################################################################

sub distribute_seq {
  my ($id, $sequence, $genome_num) = @_;
  my $placed = 0;
  my $seq_len = length $sequence;
  for (my $fidx=0; $fidx<$num_limits; $fidx++) {
    if ($seq_len <= $limits[$fidx]) {
      print {$fha[$fidx]} "$id\n$sequence";
      $placed = 1;
      last;
    }
  }
  if (!$placed) { print {$fha[$#fha]} "$id\n$sequence"; }
}
