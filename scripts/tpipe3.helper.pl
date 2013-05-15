#!/usr/bin/perl -w

my $DEBUG = 0;
my $DEBUG2 = 0;
my $MAX_NUM_LANES = 8;

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
#  November 2, 2012
#
#  tpipe3.helper.pl
#    
#     input:
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

my($raw_data_dir,$file_name_root,$work_dir,
   $dest_dir,$index_name,$sample_name,$sem);

GetOptions( 
  "raw_data_dir=s" => \$raw_data_dir,
  "file_name_root=s" => \$file_name_root,
  "work_dir=s" => \$work_dir,
  "dest_dir=s" => \$dest_dir,
  "index_name=s" => \$index_name,
  "sem=s" => \$sem,
  "sample_name=s" => \$sample_name,
);

sub printUsage{
  print "usage:\n";
  print "	-raw_data_dir : raw.data.dir.with.path  (REQUIRED)\n";
  print "	-file_name_root : filename.root.eg.C10GLACXX  (REQUIRED)\n";
  print "	-work_dir : working.dir.with.path  (REQUIRED)\n";
  print "	-dest_dir : destin.dir.for.assembled.data  (REQUIRED)\n";
  print "	-index_name : index.name  (e.g. GSLv2-7_50)";
  print "	-sem : semaphore.name  (e.g. sem.2)";
  print '	-sample_name : sample.name  (e.g. Aplysia_9thDay_veliger)';
  print "\n";  print "\n";
  print "Example call:\n\ntpipe3.helper.pl -raw_data_dir  /srv/data2/pwilliams/KH12823 -file_name_root C10GLACXX -work_dir /srv/data2/pwilliams/work_dir -dest_dir /srv/data2/pipeline -index_name GSLv2-7_50 -sample_name Aplysia_9thDay_veliger -sem sem.2 &\n";
  print "\n";
};

if(!defined($raw_data_dir)) { printUsage(); exit 0; }
if(!defined($file_name_root)) { printUsage(); exit 0; }
if(!defined($work_dir)) { printUsage(); exit 0; }
if(!defined($dest_dir)) { printUsage(); exit 0; }
if(!defined($index_name)){ printUsage(); exit 0; }
if(!defined($sem)){ printUsage(); exit 0; }
if(!defined($sample_name)){ printUsage(); exit 0; }
# print "$0 @ARGV\n";

my $log = "log." . $sample_name;
if (-e $log) { my $cmd = "rm -f $log"; system($cmd); if ( $? ) { die "Command failed: $cmd: $!"; } }
open my $fh, '>', $log or die "Unable to open $log for writing:";
# open LOG,">$log" or die "Unable to open $log for writing:";
my $ofh = select $fh; $| = 1; select $ofh; # flush ... no buffering

print $fh "raw_data_dir = $raw_data_dir\n";
print $fh "file_name_root = $file_name_root\n";
print $fh "work_dir = $work_dir\n";
print $fh "dest_dir = $dest_dir\n";
print $fh "index_name = $index_name\n";
print $fh "sample_name = $sample_name\n";
print $fh "sem = $sem\n";

print $fh "\n"; print $fh `date`; print $fh "\n";

my $cmd;

##########################################################################
# confirm all files/dirs actually exist before starting
##########################################################################

if (not -e $raw_data_dir) {
  print $fh "$raw_data_dir does not exist\n";
  exit 0;
}

if (!$DEBUG) {
    if (-e $work_dir) { 
      $cmd = "rm -f $work_dir"; 
      print $fh "$cmd\n";
      if (!$DEBUG) {
	  system($cmd);
	  if ( $? ) { die "Command failed: $cmd: $!"; }
      }
    }
    mkdir "$work_dir" or die "unable to mkdir $work_dir:";
    print $fh "mkdir $work_dir\n";

  if (not -e $dest_dir) {
    mkdir "$dest_dir" or die "unable to mkdir $dest_dir:";
    print $fh "destination dir does not exist $dest_dir, creating it\n";
    print $fh "mkdir $dest_dir\n";
    }
}

my @files_to_delete = "";

start_timing(101);

my $num_files = 0;
my @gz_data_files;


##########################################################################
# mk work_dir & change to it.
##########################################################################

print $fh "\nchdir $work_dir\n\n";
chdir $work_dir or die "Unable to chdir $work_dir:";

##########################################################################
# read raw data dir and put all files with this index_name into 
# gz_data_files array.
##########################################################################

opendir(DIR, "$raw_data_dir") or die "unable to opendir $raw_data_dir: $!";
while (my $f = readdir(DIR)) {
  next unless ($f =~ m/$index_name/);
  push @gz_data_files, $f;
#    print $fh "pushed $f onto gz_data_files array\n";
  $num_files++;
}
closedir(DIR)  or die "unable to closedir $raw_data_dir: $!";

if ($num_files == 0) { 
    print $fh "Found $num_files raw data files for $index_name; skipping this index_name.\n";
    next;
} else {
    print $fh "Found $num_files raw data files for $index_name\n";
}

##########################################################################
# parse each gz data file name to find lane num, pair num and enter file
# name into 2D array (gz_data_array), indexed by lane and pair num.
##########################################################################

my @gz_data_array;
my $last_lane_num;

# initialize 2D array
foreach my $n (1..$MAX_NUM_LANES) {
  my @tmp = "";
  $gz_data_array[$n] = [@tmp];
}

# put file name into 2D array
my $num_indiv_raw_files = 0;
foreach my $f (@gz_data_files) {
  if (!($f =~ /($file_name_root)_s(\d)_(\d)_($index_name)\.*/)) {
    print $fh "unable to parse gz file: $f to find $file_name_root\n";
    exit 0;
  }
  my $lane_num = $2;
  my $pair_num = $3;
  $gz_data_array[$lane_num][$pair_num] = $f;
  $num_indiv_raw_files++;
}
my $num_lanes_used = $num_indiv_raw_files/2;

if ($DEBUG) {
print $fh "==========================================\ngz_data_array:\n";
foreach my $n (1..$MAX_NUM_LANES) {
  foreach my $n2 (1..2) {
      if (exists $gz_data_array[$n][$n2]) {
	print $fh "$gz_data_array[$n][$n2]\n";
      }
  }
}
print $fh "==========================================\n";
}

##########################################################################
# copy gz files to work dir and gunzip them.
##########################################################################

start_timing();
foreach my $gz_file (@gz_data_files) {
    my $gz_file_with_path = $raw_data_dir . '/' . $gz_file;
    $cmd = "cp $gz_file_with_path $work_dir";
    print $fh "$cmd\n";
    if (!$DEBUG) {
	system($cmd);
	if ( $? ) { die "Command failed: $cmd: $!"; }
    }
    $cmd = "gunzip $gz_file";
    print $fh "$cmd\n";
    if (!$DEBUG) {
	system($cmd);
	if ( $? ) { die "Command failed: $cmd: $!"; }
    }
    push @files_to_delete, $gz_file;
}
end_timing($fh,"Time to gunzip files:");
print $fh "\n";

##########################################################################
# put fastq file names into 2D array.
##########################################################################

my @fq_data_array;
my @all_fq_files;

foreach my $lane_num (1..$MAX_NUM_LANES) {
  my @tmp = "";
  $fq_data_array[$lane_num] = [@tmp];
}
my $num_gz_files_pushed = 0;
foreach my $lane_num (1..$MAX_NUM_LANES) {
  foreach my $pair_num (1..2) {
      if (exists $gz_data_array[$lane_num][$pair_num]) {
	 $gz_data_array[$lane_num][$pair_num] =~ /(.*)\.gz/;
	 $fq_data_array[$lane_num][$pair_num] = $1;
	 push @all_fq_files, $1;
	 push @files_to_delete, $1;
	 $num_gz_files_pushed++;
      }
  }
}  
if($num_indiv_raw_files != $num_gz_files_pushed) {
  print $fh "$num_indiv_raw_files != $num_gz_files_pushed\n";
  exit 0;
}

##########################################################################
# print $fh names of all files in work_dir
##########################################################################
=stop
opendir(DIR, "$work_dir") or die "unable to opendir $work_dir: $!";
my @files_in_dir = readdir(DIR);
closedir(DIR)  or die "unable to closedir $work_dir: $!";
print $fh "files in $work_dir : @files_in_dir\n";
=cut
##########################################################################
# run cutadapt on all fastq files
##########################################################################

start_timing();
foreach my $in_file (@all_fq_files) {
    my $out_file = $in_file . ".cutadpated";
    $cmd = "(cutadapt -q 20 -m 1 $in_file > $out_file) >& /dev/null";
    print $fh "$cmd\n";
    if (!$DEBUG) {
	system($cmd);
	if ( $? ) { die "Command failed: $cmd: $!"; }
    }
    push @files_to_delete, $out_file;
}
end_timing($fh,"time for cutadapt:");
print $fh "\n";

##########################################################################
# interleave each pair of cutadapted fastq files
##########################################################################

my @interleaved_fq_files;

start_timing();
foreach my $lane_num (1..$MAX_NUM_LANES) {
  if (exists $fq_data_array[$lane_num][1]) {
    my $out_file = $sample_name . "_lane_" . $lane_num . "_interleaved.fastq";
    my $f1 = $fq_data_array[$lane_num][1] . ".cutadpated";
    my $f2 = $fq_data_array[$lane_num][2] . ".cutadpated";
    $cmd = "python /home/pwilliams/python_software/files/scripts/InterleaveFastq.py -f1 $f1 -f2 $f2 -o $out_file >& /dev/null";

    print $fh "$cmd\n";
    if (!$DEBUG) {
      system($cmd);
      if ( $? ) { die "Command failed: $cmd: $!"; }
     }
     push @interleaved_fq_files, $out_file;
    push @files_to_delete, $out_file;
  }
}
end_timing($fh,"time for interleave:");
print $fh "\n";

##########################################################################
# convert fastq files to fasta format
##########################################################################

#   fq_all2stdlpl needs extra arg if have \1 and \2

my @interleaved_fasta_files;
start_timing();
foreach my $file (@interleaved_fq_files) {
    $file =~ /(.*)\.fastq/;
    my $outfile = $1 . ".fasta";
    $cmd = "(perl /srv/data2/software/PerlScripts/Format/fq_all2std.pl fq2fa $file > $outfile) >& /dev/null";
    print $fh "$cmd\n";
    if (!$DEBUG) {  
	system($cmd);
	if ( $? ) { die "Command failed: $cmd: $!"; }
    }
    push @interleaved_fasta_files, $outfile;
    push @files_to_delete, $outfile;
}
end_timing($fh,"time for conversion from fastq to fasta format:");
print $fh "\n";

##########################################################################
# print $fh names of all files in work_dir
##########################################################################
=stop
opendir(DIR, "$work_dir") or die "unable to opendir $work_dir: $!";
@files_in_dir = readdir(DIR);
closedir(DIR)  or die "unable to closedir $work_dir: $!";
print $fh "files in $work_dir : @files_in_dir\n";
=cut

##########################################################################
# cat together interleaved data
##########################################################################

start_timing();
my $combined_interleaved_fasta_data = $sample_name . "_combined_interleaved.fasta";
$cmd = "touch $combined_interleaved_fasta_data";
print $fh "$cmd\n";
if (!$DEBUG) {  
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

foreach my $file (@interleaved_fasta_files) {
    $cmd = "cat $file >> $combined_interleaved_fasta_data";
    print $fh "$cmd\n";
    if (!$DEBUG) {
	system($cmd);
	if ( $? ) { die "Command failed: $cmd: $!"; }
    }
}
push @files_to_delete, $combined_interleaved_fasta_data;
end_timing($fh,"time to cat all files together:");
print $fh "\n";

##########################################################################
# run khmer ....
##########################################################################

# scp -P 32080 Aplysia_9thDay_veliger_interleaved.fasta administrator@128.227.123.35:~

start_timing();
my $khmer_out = $combined_interleaved_fasta_data . ".keep";
$cmd = "/srv/data2/software/Khmer/scripts/normalize-by-median.py -C 30 -k 20 -N 1 -x 4.0e9 $combined_interleaved_fasta_data  >& /dev/null";


print $fh "$cmd\n";
if (!$DEBUG) {
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
}
push @files_to_delete, $khmer_out;
end_timing($fh,"time to run khmer:");
print $fh "\n";

##########################################################################
# deinterleave data
##########################################################################

my $deinterleaved1 = $khmer_out . ".1";
my $deinterleaved2 = $khmer_out . ".2";

start_timing();
$cmd = "python /home/pwilliams/python_software/files/scripts/seqtools.py -f $khmer_out -ft fasta --split-paired-end >& /dev/null";
print $fh "$cmd\n";
if (!$DEBUG) {
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
}
end_timing($fh,"time to run deinterleave:");
push @files_to_delete, $deinterleaved1;
print $fh "\n";

##########################################################################
# run trinity assembler
##########################################################################

start_timing();
my $trinity_out_dir = "trinity_out";
my $trinity_out = $trinity_out_dir . "/Trinity.fasta";
$cmd = "perl /srv/data2/software/Trinity/Trinity.pl --seqType fa --left $deinterleaved1 --right $deinterleaved2 --JM 100G --inchworm_cpu 20 --CPU 20 --output $trinity_out_dir > /dev/null";
print $fh "$cmd\n";
if (!$DEBUG) {
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
}
end_timing($fh,"time to run Trinity on $sample_name:");
print $fh "\n";

##########################################################################
# move assembled data to blast pipe location
##########################################################################

$cmd = "ls -lh $trinity_out_dir";
print $fh "$cmd\n";
#  if (!$DEBUG) {
  open (GR, "$cmd |") or die "unable to open  pipe:";
  while (my $line = <GR>) { print $fh "$line"; } close GR;
#  }
print $fh "\n";

my $tdir = $dest_dir . "/" . $sample_name;
$cmd = "mkdir $tdir";
print $fh "$cmd\n";
if (!$DEBUG) {
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
}

my $tgt = $dest_dir . "/" . $sample_name . "/" . $sample_name . ".fasta";
$cmd = "mv $trinity_out $tgt";
print $fh "$cmd\n";
if (!$DEBUG) {
    system($cmd);
    if ( $? ) { die "Command failed: $cmd: $!"; }
}

print $fh "\n";

##########################################################################
# remove work dir
##########################################################################

opendir(DIR, "$work_dir") or die "unable to opendir $work_dir: $!";
my @files_in_dir = readdir(DIR);
closedir(DIR)  or die "unable to closedir $work_dir: $!";
print $fh "files in $work_dir : @files_in_dir\n\n";

$cmd = "rm -fr $trinity_out_dir";
print $fh "$cmd\n";
if (!$DEBUG) {
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

print $fh "\nchdir ..\n";
chdir ".." or die "Unable to chdir ..:";

$cmd = "rm -fr $work_dir";
print $fh "$cmd\n";
if (!$DEBUG) {
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

print $fh "\n";

print $fh "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
end_timing($fh,"Total time to run pipeline for $index_name -- $sample_name :",101);
print $fh "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";

$cmd = "touch $sem";
print $fh "$cmd\n";
if (!$DEBUG) {
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

my @start_time;

sub start_timing {
  my ($N) = @_;
  if (not defined $N) { $N = 0; }
  $start_time[$N] = time;
#  print $fh_date();
}

sub end_timing {
  my ($fh,$msg, $N) = @_;
  if (not defined $N) { $N = 0; }
  my $etime = ((time - $start_time[$N]) / 60.0);
#  print $fh_date();
  my $dat = `date`;
  chomp $dat;
  printf $fh "$msg:\t%8.3f\tmins\t(\t%8.3f\thrs)\t$dat\n", $etime, $etime/60.0;
  return $etime;
}

