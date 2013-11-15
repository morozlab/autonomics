#!/usr/bin/perl -w
use strict;
use Sys::Hostname;
use File::Basename;
use Cwd;
use Getopt::Long;

my ($query, $arg1) = @ARGV;
if ((not defined $query)) {
    print "\nUsage: $0 <cmd pa | pn | pname | pid | jobs | q | run |nf | nfns | snf >  <arg1 optional>\n";
    print " 	pipe_status pa: prints <project_id> <project_name> for all projects in database, finished and unfinished, sorted alphabetically\n";
    print " 	pipe_status pn: prints <project_id> <project_name> for all projects in database, finished and unfinished, sorted by project_id\n";
    print " 	pipe_status pname <proj_id>: prints <project_name for that <project_id>\n";
    print " 	pipe_status pid <proj_name>: prints <project_id for that <project_name>\n";
    print " 	pipe_status jobs <proj_id>: prints status of all jobs for that <project_id>\n";
    print " 	pipe_status q : prints jobs in queue\n";
    print " 	pipe_status run : prints all jobs that are currently running\n";
    print " 	pipe_status nf : prints all jobs that are not finished\n";
    print " 	pipe_status snf : prints all jobs that are started but not finished\n";
    print " 	pipe_status nfns : prints all jobs that are not finished and not yet started\n";
    exit 0;
}

my $end =  "\ | mysql -u zeroclick --password=Whitney2011 --database=zero_click -vvv -t";
my $tend =  " | mysql -u zeroclick --password=Whitney2011 --database=zero_click";

if ($query eq 'pid') {
    if (not defined $arg1) {
	print "need arg1 = project_name\n";
	exit 0;
    }
  my $cmd = "echo \"select project_id,project_name from pn_mapping where project_name = \'$arg1\'\"" . $tend;
  # print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'pname') {
    if (not defined $arg1) {
	print "need arg1 = pid\n";
	exit 0;
    }
  my $cmd = "echo \"select project_id,project_name from pn_mapping where project_id = \'$arg1\'\"" . $tend;
#  print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if (($query eq 'pa')) {
  my $cmd = "echo select project_id,project_name from pn_mapping" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'pn') {
  my $cmd = "echo select project_id,project_name from  pn_mapping order by project_id" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'nf') {
  my $cmd = "echo \"select * from jn_mapping where finished = \'N\'\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'nfns') {
  my $cmd = "echo \"select * from jn_mapping where started = \'N\' and finished = \'N\'\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'snf') {
  my $cmd = "echo \"select * from jn_mapping where started = \'Y\' and finished = \'N\'\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'run') {
  my $cmd = "echo \"select * from jn_mapping where started = \'Y\' and finished = \'N\'\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'err') {
  my $cmd = "echo \"select * from jn_mapping where finished = \'E\'\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'q') {
  my $cmd = "echo \"select * from quenew\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'jobs') {
  if (not defined $arg1) {
      print "need arg1 = pid\n";
      exit 0;
  }
  my $cmd = "echo \"select * from jn_mapping where project_id = $arg1\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}

if ($query eq 'star') {
  if (not defined $arg1) {
      print "need arg1 = table_name\n";
      exit 0;
  }
  my $cmd = "echo \"select * from $arg1\"" . $tend;
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}


# select job_id from args where pipeline_args ='--assembler mira';




if ($query eq 'dt') {
    if (not defined $arg1) {
	print "need arg1 = table_name\n";
	exit 0;
    }
  my $cmd = "echo \"describe $arg1\"" . $end;
  print "$cmd\n";
  system($cmd);
  if ( $? ) { die "Command failed: $cmd: $!"; }
}
