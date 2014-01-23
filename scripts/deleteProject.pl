#!/usr/bin/perl -w
use warnings;
use strict;
use DBI;
use File::Basename;
use File::Path;
# use Cwd;
use POSIX;

my $apath  = $ENV{AUTONOMICS_PATH };
my $config_file = $apath . '/config.pl';
my %config = do $config_file;
my $root_pw = $config{db_root_password};

my $db_dir = $ENV{NEUROBASE_SEQ_PATH};

my ($pid) = @ARGV;
if (not defined $pid) {
    print "\nUsage: $0 <project_id>\n";
    exit 0;
}

if (not isdigit $pid) { 
   print "projectID: $pid  must be an integer \n"; 
   exit 0;
}

my $dsn = "dbi:mysql:database=moroz_lab;host=localhost";
my $dbh = DBI->connect($dsn, "root", $root_pw);

print "Deleting database files:  fasta file and formatdb files\n";
my $cmd = "/bin/rm -fr $db_dir/$pid";
print "$cmd\n";
system ($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

print "Removing annotation alignments\n";
my $query = "SELECT * FROM $pid\_sequences";
my $sth = $dbh->prepare($query);
$sth->execute();
my $seqTableExists = 0;
while(my $row = $sth->fetchrow_hashref()){
  my $sbID = $row->{sb_id};
  my $removal = "DELETE FROM annotation_alignments WHERE sb_id='" . $sbID . "'";
  $dbh->do($removal);
  $seqTableExists = 1;
}

my $project_selector = " WHERE project_id='$pid'";

#remove entries from the project_files table  project_files not used
print "Removing project files entries\n";
$query = "DELETE FROM project_files WHERE pid ='$pid'";
$dbh->do($query);

$query = "DELETE FROM load_info WHERE project_id ='$pid'";
$dbh->do($query);

print "Removing sb_catalog entries\n";
#remove entries from the sb_catalog table
$query = "DELETE FROM sb_catalog WHERE projectID ='$pid'";
$dbh->do($query);

print "Removing homology entries\n";
#remove entries from the homology tables
$query = "DELETE FROM best_annotations" . $project_selector;
$dbh->do($query);
$query = "DELETE FROM homology" . $project_selector;
$dbh->do($query);
$query = "DELETE FROM sorted_homology" . $project_selector;
$dbh->do($query);

print "Removing GO_annotation entries\n";
#remove entries from the GO annotation table
$query = "DELETE FROM go_annotation" . $project_selector;
$dbh->do($query);

#remove entries from the new GO annotation table
$query = "SELECT * FROM go_annotation_new, $pid\_sequences WHERE go_annotation_new.sb_id = $pid\_sequences.sb_id LIMIT 0,1";
$sth = $dbh->prepare($query);
$sth->execute();
my $has_new_go = 0;
while(my $ref = $sth->fetchrow_hashref()){
  $has_new_go = 1;
}
if($has_new_go){
  print "This project has new style GO-annotation storage, removing.";
  $query = "DELETE FROM go_annotation_new USING go_annotation_new, $pid\_sequences where go_annotation_new.sb_id = $pid\_sequences.sb_id";
  $dbh->do($query);
}

print "Removing cross_comparison entries\n";
#remove entries from the cross_comparisons table
$query = "DELETE FROM cross_comparisons WHERE project_id_1='$pid' OR project_id_2='$pid'";
$dbh->do($query);

print "Removing cross_comparisons_list entries\n";
#remove entries from the cross_comparisons_list table
$query = "DELETE FROM cross_comparisons_list WHERE project_id_1='$pid' OR project_id_2='$pid'";
$dbh->do($query);

print "Removing go_categories entries\n";
#remove entries from the go_categories table
$query = "DELETE FROM go_categories" . $project_selector;
$dbh->do($query);

print "Reomving kegg entries.\n";
$query = "DELETE FROM kegg_annotations" . $project_selector;
$dbh->do($query);

print "Removing pfam entries.\n";
$query = "DELETE FROM pfam_domain_counts" . $project_selector;
$dbh->do($query);
$query = "DELETE FROM pfam_annotations". $project_selector;
$dbh->do($query);

print "Emptying sequence table\n";
#empty the sequence table for this project
#if($seqTableExists){
#  $query = "TRUNCATE TABLE $pid\_sequences";
#  $dbh->do($query);
  #delete the sequence table for this project
  $query = "DROP TABLE IF EXISTS $pid\_sequences";
  $dbh->do($query);
# }

#remove this project from the project directory table

$query = "DELETE FROM project_directory WHERE projectID='$pid'";
$dbh->do($query);

print "DELETE COMPLETED\n";

