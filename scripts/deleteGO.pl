#!/usr/bin/perl -w
use warnings;
use strict;
use DBI;
use Sys::Hostname;
use File::Basename;
use File::Path;
use Cwd;

my ($projectID) = @ARGV;
if ((not defined $projectID)) {
    print "\nUsage: $0 <project_id>\n";
    exit 0;
}

my $host = $ENV{HOST};
if ($host =~ /oem/) {  $host = "oem"; }
if ($host =~ /hpc/) {  $host = "hpc"; }
if ($host =~ /acis/) {  $host = "acis"; }

my $dsn = "dbi:mysql:database=moroz_lab;host=localhost";
my $dbh = DBI->connect($dsn, "root", "Moof2011");

my $db_dir = "";
if ($host eq "hpc") {
    $db_dir = "/var/www/html/neurobase/seq_view/database/";
}
if ($host eq "oem") {
    $db_dir = "/var/www/seq_view/database";
}
if ($host eq "acis") {
    $db_dir = "/home/pwilliams/foo";
}

my $project_selector = " WHERE project_id='$projectID'";

print "removing entries from the new GO annotation table\n";
my $query = "SELECT * FROM go_annotation_new, $projectID\_sequences WHERE go_annotation_new.sb_id = $projectID\_sequences.sb_id LIMIT 0,1";
my $sth = $dbh->prepare($query);
$sth->execute();
my $has_new_go = 0;
while(my $ref = $sth->fetchrow_hashref()){
  $has_new_go = 1;
}
my $numd = 0;
if($has_new_go){
    $numd++;
  $query = "DELETE FROM go_annotation_new USING go_annotation_new, $projectID\_sequences where go_annotation_new.sb_id = $projectID\_sequences.sb_id";
  $dbh->do($query);
}
print "Removing go_categories entries\n";
#remove entries from the go_categories table
$query = "DELETE FROM go_categories" . $project_selector;
$dbh->do($query);
print "done\n";
