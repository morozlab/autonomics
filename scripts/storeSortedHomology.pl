use warnings;
use strict;
use DBI;
# use File::copy;

my $apath  = $ENV{AUTONOMICS_PATH } . "/";
my $scriptpath = $apath . "scripts/";

my $tmpPath  = $ENV{NEUROBASE_LOAD_DATA_PATH } . "/";

my ($projName, $projectID, $sort, $sqldb, $user, $passwd, $debug) = @ARGV;
if ((not defined $projectID) ||
    (not defined $sort) ||
    (not defined $sqldb) ||
    (not defined $user) ||
    (not defined $passwd) ||
    (not defined $projName) ||
    (not defined $debug)) 
{
    print "\nUsage: $0 <proj_name> <proj_id> <sort> <user> <passwd> <debug>\n";
    print "  <sort>  type: 1 for evalue, 2 for abundance>\n";
    exit 0;
}

print "$0 $projName $projectID $sort $debug\n";

print "dbi:mysql:database=".$sqldb.";host=localhost\n";
#my $dsn = "dbi:mysql:database=moroz_lab;host=localhost";

my $dsn = "dbi:mysql:database=".$sqldb.";host=localhost";
my $dbh = DBI->connect($dsn, $user, $passwd);

$tmpPath = $tmpPath . '/' . $projName . "/";;
print "tmpPath: $tmpPath\n";
my $outfile = $tmpPath . "homologysort.txt";
print "outfile = $outfile\n";

# out1 & out2 only used if debug = 1
my $out1 = $scriptpath . "homologysort1.txt" . $$;
my $out2 = $scriptpath . "homologysort2.txt" . $$;



open(OUTFILE, ">$outfile") or die "Could not open temp output file: $outfile!\n";

print " opened $outfile for writing\n";

if($sort == 1){
    if ($debug) {  open(OUT1, ">$out1") or die "Could not open temp output file: $out1\n"; }
    if ($debug) {  print "storeSortedHomology.pl: sorting homology table by evalue into: $outfile\n";}
    if ($debug) {  print "cols of $outfile are:\n\tpid, 1, ranking, sb_id, annot_id, evalue, abund, source\n";}
  my $query = "SELECT homology.sb_id, annotation_id, evalue, abundance, source FROM homology, $projectID\_sequences WHERE homology.project_id='" . $projectID . "' AND homology.sb_id = $projectID\_sequences.sb_id ORDER BY evalue";
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my $ranking = 0;
  while(my $row = $sth->fetchrow_hashref()){
    if($ranking == 0){
      print OUTFILE $projectID . "\t" . "1\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";
      if ($debug) {      print OUT1 $projectID . "\t" . "1\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t"; }
    }
    else{
      print OUTFILE "\n";
      print OUTFILE $projectID . "\t" . "1\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";
      if ($debug) {      print OUT1 "\n"; }
      if ($debug) {      print OUT1 $projectID . "\t" . "1\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";}
    }
    $ranking++;
  }
    if ($debug) {  close(OUT1); }

}
elsif($sort == 2){
  #can only run after the project has been sorted by evalue
    if ($debug) {  open(OUT2, ">$out2") or die "Could not open temp output file: $out2\n"; }
    if ($debug) {  print "storeSortedHomology.pl: sorting all seqs in homology table with abundance > 1 by abundance into $outfile\n";}
    if ($debug) {  print "cols of $outfile are:\n\tpid, 2, ranking, sb_id, annot_id, evalue, abund, source\n";}

  my $query = "SELECT homology.sb_id, annotation_id, evalue, $projectID\_sequences.abundance, source FROM homology, " . $projectID . "_sequences WHERE homology.sb_id = " . $projectID . "_sequences.sb_id AND homology.project_id ='" . $projectID . "' AND " . $projectID . "_sequences.abundance > 1 ORDER BY " . $projectID . "_sequences.abundance DESC";
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my $ranking = 0;
  while(my $row = $sth->fetchrow_hashref()){
    if($ranking == 0){
      print OUTFILE $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";
      if ($debug) {      print OUT2 $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";}
    }
    else{
      print OUTFILE "\n";
      print OUTFILE $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";
      if ($debug) {      print OUT2 "\n";}
      if ($debug) {      print OUT2 $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t" . $row->{abundance} . "\t" . $row->{source} . "\t";}
    }
    $ranking++;
  }

    if ($debug) {  print "storeSortedHomology.pl: sorting all seqs in sorted_homology table with sort_id = 1 & abundance = 1 by ranking & appending them to $outfile\n";}
    if ($debug) {  print "those that are appended have 1 in abuindance col\n";}

    if ($debug) {   print OUT2 "==============================\n";}

  $query = "SELECT sorted_homology.sb_id, annotation_id, evalue, " . $projectID . "_sequences.abundance, source FROM sorted_homology, " . $projectID . "_sequences WHERE sorted_homology.sb_id = " . $projectID . "_sequences.sb_id AND sorted_homology.project_id ='" . $projectID . "' AND sort_id='1' AND " . $projectID . "_sequences.abundance = '1' ORDER BY ranking";

  $sth = $dbh->prepare($query);
  $sth->execute();
  while(my $row = $sth->fetchrow_hashref()){
    print OUTFILE "\n";
    print OUTFILE $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t1\t" . $row->{source} . "\t";
    if ($debug) {    print OUT2 "\n";}
    if ($debug) {    print OUT2 $projectID . "\t" . "2\t" . "$ranking\t" . $row->{sb_id} . "\t" . $row->{annotation_id} . "\t" . $row->{evalue} . "\t1\t" . $row->{source} . "\t";}
    $ranking++;
  }
    if ($debug) {  print "done.\n";}
  #update the project_directory to say that this project has abundance data
  $query = "UPDATE project_directory SET has_abundance ='Y' WHERE projectID='$projectID'";
  $dbh->do($query);
    if ($debug) {  close(OUT2);}
}

close(OUTFILE);

# exit 0;

#load the data in homology file

my $query = "";

# if ($debug) {  print 'LOAD DATA INFILE ' . $dbh->quote($outfile) . ' REPLACE INTO TABLE sorted_homology';}

my $cmd = "python load.homology.py --user $user --password $passwd --outfile $outfile";
print "$cmd\n";
system($cmd);
if ( $? ) { die "Command failed: $cmd: $!"; }

=stop
my $result = eval {
   $query = "LOAD DATA INFILE " . $dbh->quote($outfile) . " REPLACE INTO TABLE sorted_homology";
   print "QUERY: $query\n";
   $dbh->do("LOAD DATA INFILE " . $dbh->quote($outfile) . " REPLACE INTO TABLE sorted_homology");
   print "success\n";
};
unless ($result) {
   $query = "LOAD DATA LOCAL INFILE " . $dbh->quote($outfile) . " REPLACE INTO TABLE sorted_homology";
   print "Q2: $query\n";
   $dbh->do("LOAD DATA LOCAL INFILE " . $dbh->quote($outfile) . " REPLACE INTO TABLE sorted_homology");
   print "success\n";
}
=cut

unlink($outfile);

my $t1 = $scriptpath . "ss.done1";
my $t2 = $scriptpath . "ss.done2";

print "t1: $t1\n";
print "t2: $t2\n";

if($sort == 1) {
  my $cmd = "touch $t1";
  if ($debug) {  print "$cmd\n";}
  system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }
} else {
  my $cmd = "touch $t2";
  if ($debug) {  print "$cmd\n";}
  system($cmd);
 if ( $? ) { die "Command failed: $cmd: $!"; }
}

print "storeSortedHomology() done for sort# $sort\n";
