#!/usr/bin/perl -w

my $DEBUG = 0;
my $DEBUG2 = 0;
my $MAX_NUM_LANES = 8;
my $MAX_NUM_PROCS = 3;

use strict;
use IO::File;
use File::Basename;
use Time::HiRes qw( gettimeofday usleep tv_interval );
use Time::localtime;
use Sys::Hostname;
use File::Path;
use Getopt::Long;

##############################################################################
#
#  Peter L. Williams
#  October 26, 2012
#
#  tpipe3.pl
#
#     PREPROCESS THE SPREADSHEET FILE AS FOLLOWS:
#       - save spreadsheet from excel as TAB DELIMITED file, call it FTAB
#       - run " tr '\r' '\n' FTAB > spread
#       - edit spread & remove first & last lines & any other non useful lines.
#       - also check all column 4 names to be sure they are meaningful --- 
#       -   'unidentified species' is not meaningful - col 4 entry used as
#       -   name for final assembled data.  To enter a tab in emacs ^q<tab>.
#
#     COPY SPREADSHEET AND RAW DATA TO ACIS:
#       - cd /media/My\ Passport/
#       - scp -r ./KH12823 pwilliams@128.227.70.246:~pwilliams/data2
#       - cd dir with modified spreadsheet 
#       - scp spreadsheet pwilliams@128.227.70.246:~pwilliams/data2/KH12823
#
#
#     scp raw data to acis:  
#     input: - raw data dir (dir, with path, that holds data from hard disk)
#            - readme file with path (normally in above dir)
#            - raw data file name root (e.g. 'C10GLACXX')
#            - spreadsheet with path (from Andrea) in tab delim format with ^M
#              removed: to remove ^M: 
#              " tr '\r' '\n' < spreadshhet.prelim.tab > spreadsheet.tab "
#            - path where can create working dir for intermediate data.
#            - target dir for assembly output.
#            - start_index_name { e.g. GSLv2-7_50 } [optional: 
#              use to restart pipeline]
#
#     output: assembled data put in raw_data_dir in subdir named <sample name>
#             (spreadsheet col 4).
#         
#     runs:  - reads raw data dir and for each data set in that dir:
#              -- determines how many lanes used for sequencing that data set
#              -- copies raw_data set to work dir
#              -- gunzip on raw data in work dir
#              -- cutadapt -q 20 -m 1
#              -- interleaves each pair of fastq files (InterleaveFastq.py)
#              -- converts interleaved fastq files to fasta format
#                 (fq_all2std.pl fq2fa)
#              -- concatenates all fasta files into one file.
#              -- khmer to normalize concatenated data 
#                 (normalize-by-median.py -C 30 -k 20 -N 1 -x 4.0e9)
#              -- deinterleaves above result 
#                 (SequenceUtils.py -f outfile -ft fasta  --split-paired-end)
#              -- assembler: (
#                 (Trinity.pl  --seqType fa --left <deinterleaved1.data> 
#                         --right <deinterleaved2.data> --JM 100G 
#                         --inchworm_cpu 20 --CPU 20 --output <trinity_out>)
#
#              -- copies trinity output to target dir & to blast pipeline dir.
#              -- removes files in work_dir
#
#           - each step is timed and log info printed for each step.
#
#           - run with output to log, i.e. >& log.run &
#            
##############################################################################

my($raw_data_dir,$readme,$spreadsheet,$file_name_root,$work_dir_path,
   $dest_dir,$start_index_name,$num_procs);

GetOptions( 
  "raw_data_dir=s" => \$raw_data_dir,
  "readme=s" => \$readme,
  "spreadsheet=s" => \$spreadsheet,
  "file_name_root=s" => \$file_name_root,
  "work_dir_path=s" => \$work_dir_path,
  "dest_dir=s" => \$dest_dir,
  "num_procs=i" => \$num_procs,
  "start_index_name=s" => \$start_index_name,
);

sub printUsage{
  print "usage:\n";
  print "	-raw_data_dir : raw.data.dir.with.path  (REQUIRED)\n";
  print "	-readme : readme.file.with.path  (REQUIRED)\n";
  print "	-spreadsheet : spread.sheet.with.path.tab.delim  (REQUIRED)\n";
  print "	-file_name_root : filename.root.eg.C10GLACXX  (REQUIRED)\n";
  print "	-work_dir_path : path to where we can create a working_dir  (REQUIRED)\n";
  print "	-dest_dir : destin.dir.for.assembled.data  (REQUIRED)\n";
  print "	-num_procs : numb jobs to run concurrently (OPTIONAL,  DEFAULT: 3 - max is $MAX_NUM_PROCS)\n";
  print '	-start_index_name : start.index.name  (OPTIONAL, e.g. GSLv2-7_50   DEFAULT: "")';
  print "\n";  print "\n";
  print "Example call:\n\ntpipe3.pl -raw_data_dir  /srv/data2/pwilliams/KH12823 -readme /srv/data2/pwilliams/KH12823/README_C10GLACXX -spreadsheet /srv/data2/pwilliams/KH12823/spread.tabs -file_name_root C10GLACXX -work_dir_path /srv/data2/pwilliams -dest_dir /srv/data2/pipeline >& log.tpipe.12 &\n";
  print "\n";
};

if(!defined($raw_data_dir)) { printUsage(); exit 0; }
if(!defined($readme)) { printUsage(); exit 0; }
if(!defined($spreadsheet)) { printUsage(); exit 0; }
if(!defined($file_name_root)) { printUsage(); exit 0; }
if(!defined($work_dir_path)) { printUsage(); exit 0; }
if(!defined($dest_dir)) { printUsage(); exit 0; }
if(!defined($num_procs)) { $num_procs = 3; }
if(!defined($start_index_name)) { $start_index_name = ""; }
# print "$0 @ARGV\n";

print "raw_data_dir = $raw_data_dir\n";
print "readme = $readme\n";
print "spreadsheet = $spreadsheet\n";
print "file_name_root = $file_name_root\n";
print "work_dir_path = $work_dir_path\n";
print "dest_dir = $dest_dir\n";
print "num_procs = $num_procs\n";
print "start_index_name = $start_index_name\n";
print "\n"; print `date`; print "\n";

my $cmd;

##########################################################################
# confirm all files/dirs actually exist before starting
##########################################################################

if (($num_procs <= 0) || ($num_procs > $MAX_NUM_PROCS)) {
    print "$num_procs must > 0 and <= $MAX_NUM_PROCS\n";
    exit 0;
}

if (not -e $raw_data_dir) {
  print "$raw_data_dir does not exist\n";
  exit 0;
}
if (not -e $readme) {
  print "$readme does not exist\n";
  exit 0;
}
if (not -e $spreadsheet) {
  print "$spreadsheet does not exist\n";
  exit 0;
}


if ((not -e $work_dir_path) && (not -d $work_dir_path)) {
  print "$work_dir_path does not exist or is not a dir\n";
  exit 0;
}  

if (!$DEBUG) {
  if (not -e $dest_dir) {
    mkdir "$dest_dir" or die "unable to mkdir $dest_dir:";
    print "destination dir does not exist $dest_dir, creating it\n";
    print "mkdir $dest_dir\n";
    }
}

##########################################################################
# read README file and enter all index names into index_names array and 
# create dictionary: machine_num (eg 1823-KH-001) to index_name. GSLv2-7_50
##########################################################################

my @index_names_with_dupes;
my $num_index_names_with_dupes = 0;
my %index_name_to_machine_num;
# my %machine_num_to_index_name;
open README, "$readme" or die "Unable to open $readme:";
while (my $line = <README>) {
  chomp $line;
# PS8704 contains:
#   Index       Sequence   Library	Name
# GSLv2-7_49    TGGGAGT    SL19946	1823-KH-0001
#  col1          col2       col3          col4
  my ($col1, $col2, $col3, $col4) = split( /\s+/, $line );
  next if (($line =~ /^\s*$/) or ($col1 eq "Index") or (!defined $col4));
  push @index_names_with_dupes, $col1;
  $num_index_names_with_dupes++;
  $index_name_to_machine_num{$col1} = $col4;
  if ($DEBUG2) {  print "$col1 ==> $col4\n"; }
#  $machine_num_to_index_name{$col4} = $col1;
}
close README;

##########################################################################
# Remove duplicate index names from list of index names and print data
##########################################################################

my %temp = ();
my @index_names = grep ++$temp{$_} < 2, @index_names_with_dupes;
my $num_index_names = @index_names;
my $num_dupes = $num_index_names_with_dupes - $num_index_names;
print "Found $num_index_names_with_dupes index names in $readme\n";
print "of which $num_dupes were duplicates, so $num_index_names unique\n";


foreach my $f (@index_names) {
  if ($DEBUG2) { print "$f\n"; }
}

##########################################################################
# Create dictionary of machine_num (eg 1823-KH-001) to 
# sample name (eg Aplysia_3d day-3)
##########################################################################

my %machine_num_to_sample_name;
my $num_sample_names = 0;
my @sample_names;
open SPREADSHEET, "$spreadsheet" or die "Unable to open $spreadsheet:";
# my $line1 = 1;
while (my $line = <SPREADSHEET>) {
#  if ($line1) { $line1 = 0; next; }
  chomp $line;
#  col1        col2      col3                     clo4
#  Illumina     #       Organism                Sample name           ...
# 1823-KH-0001	1	Aplysia california	Aplysia_4th Day early

  my ($col1, $col2, $col3, $col4, @rest) = split( /\t/, $line );
  next if (($line =~ /^\s*$/) or ($col1 eq "Index") or ($col4 eq ""));
  if (!defined $col1) { next; }
  $col4 =~ s/\s/_/g;
  $col4 =~ s/\//_/g;
  push @sample_names, $col4;
  $num_sample_names++;
  $machine_num_to_sample_name{$col1} = $col4;
  if ($DEBUG2) { print "$col1 ==> $col4\n"; }
}
close SPREADSHEET;

print "Found $num_sample_names sample names in $spreadsheet\n";
foreach my $f (@sample_names) {
    if ($DEBUG2) { print "$f\n"; }
}

##########################################################################
# initialize semaphore stuff 
##########################################################################

my @sems;
foreach my $sem (1 .. $num_procs) {
  $sems[$sem] = $work_dir_path . "/sem" . $sem;
  if (-e $sems[$sem]) {
    $cmd = "rm $sems[$sem]";
    print "$cmd\n";
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
  }
}

my @sem_used;
foreach my $i (1 .. $num_procs) {
    $sem_used[$i] = 0;
}

##########################################################################
# For each index_name: 
#     Read raw data dir to find relevant data file names and put them
#     into 2D array indexed by lane number and pair number 
#     (gz_data_array[lane_num][pair_num].  
#     Then run the pipeline on that data set.
##########################################################################

my $need_to_skip;
if ($start_index_name eq "") { $need_to_skip = 0; }
else { $need_to_skip = 1; }

my $num_data_sets_run = 0;

start_timing(201);

foreach my $index_name (@index_names) {

  print "Working on index_name = $index_name\n";

  if ($need_to_skip) {
    if ($index_name eq $start_index_name) {
      $need_to_skip = 0;
    } else {
      print "skipping $index_name because neq specified start_index $start_index_name\n";
      next;
    }
  }
   
  start_timing(101);

  ##########################################################################
  # find sample name, e.g. Aplysia_9thDay_veliger
  ##########################################################################

  my $mach_num = $index_name_to_machine_num{$index_name};
#  print "XXXX mach_num = $mach_num\n";

  if ((!defined $mach_num) or ($mach_num eq "")) { 
    next; 
    print "skipping index_name: $index_name since not on spreadsheet \n";
    # some index_names on disk do not relate to current spreadsheet
  }

  my $sample_name = 
    $machine_num_to_sample_name{$index_name_to_machine_num{$index_name}};
#  print "XXXX sample_name = $sample_name\n";
#  print "XXXX index_name = $index_name\n";

  if ((not defined $sample_name) or ($sample_name eq "")) {
      print "sample name not defined for $index_name ... skipping this index_name\n";
      next;
  }
  print "working on sample_name $sample_name\n";
  print "----------------------------------\n";

  my $work_dir = $work_dir_path ."/".$sample_name . ".work.dir";

  my $sem;
LABEL:
  my $free_sem = 0;
  foreach my $indx (1..$num_procs) {
    if (!$sem_used[$indx]) {
	$free_sem = $indx;
        last;
    }
  }
  print "free_sem = $free_sem\n";

  if ($free_sem) {   # run job
    $sem = $sems[$free_sem];
    print "sem = $sem\n";
    $sem_used[$free_sem] = 1;
    print "processing $sample_name\n";
    my $pid = fork();
    if (not defined $pid) { die "fork failed:"; }
    if ($pid == 0) {
      $cmd = "tpipe3.helper.pl -raw_data_dir  $raw_data_dir -file_name_root $file_name_root -work_dir $work_dir -dest_dir $dest_dir -index_name $index_name -sample_name $sample_name -sem $sem";
      system($cmd);
      if ( $? ) { die "Command failed: $cmd: $!"; }
      exit(0);
    } else { #parent does this
      next;  #do next index_name
    }
  } else {  # wait for sem
     my $found_a_sem = 0;
     while (1) {
       for my $i (1 .. $num_procs) {
	 if (-e $sems[$i]) {
	     $found_a_sem = 1;
	     $cmd = "rm $sems[$i]";
	     print "$cmd\n";
	     system($cmd);
	     if ( $? ) { die "Command failed: $cmd: $!"; }
	     $sem_used[$i] = 0;
	     last;
	 }
       }
       if ($found_a_sem) {
	 last;
       } else {
	 print "sleeping 1200\n";
	 sleep 1200;
       }
     }
  } #end else wait for sem
  goto LABEL;
  $num_data_sets_run++;
}

##########################################################################
# end foreach index_name loop
##########################################################################

print "\n=========================================================\n";
end_timing("Total time to run pipeline for $num_data_sets_run runs",201);
print "\n=========================================================\n\n";


my @start_time;

sub start_timing {
  my ($N) = @_;
  if (not defined $N) { $N = 0; }
  $start_time[$N] = time;
#  print_date();
}

sub end_timing {
  my ($msg, $N) = @_;
  if (not defined $N) { $N = 0; }
  my $etime = ((time - $start_time[$N]) / 60.0);
#  print_date();
  my $dat = `date`;
  chomp $dat;
  printf "$msg:\t%8.3f\tmins\t(\t%8.3f\thrs)\t$dat\n", $etime, $etime/60.0;
  return $etime;
}

sub print_date {
  open (DATE ,"date +%c |") or 
	    die "Unable to execute \"date command |\": $!";
  my $date = <DATE>;
  chomp $date;
  close (DATE) || die "unable to close DATE $!";
  print "date: $date\n";
}

