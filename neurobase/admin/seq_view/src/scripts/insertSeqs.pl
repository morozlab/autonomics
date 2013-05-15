use DBI;
use Bio::SeqIO;
use strict;
use warnings;

my $dsn = "dbi:mysql:database=moroz_lab;host=localhost";
my $dbh = DBI->connect($dsn, "root", "meow12");

my $in = shift or die;
my $tablename = shift or die;
my $fileID = shift or die;
my $type = shift or die;

print $tablename . "\n";

my $dbPath = "C:\\wamp\\www\\seq_view\\database";
#create the sequence table for the newly added sequences
my $tableStr = "show tables like '" . $tablename . "_sequences'";
my $sth = $dbh->prepare($tableStr);
$sth->execute();
my $found = 0;
while(my $hashref = $sth->fetchrow_hashref()){
  $found = 1;
  last;
}

my $fullname = $tablename . "_sequences";

#create the table to hold the sequences
if($found == 0){
  $tableStr = "CREATE TABLE " . $fullname . " (seq_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY, sb_id INT(20), project_id INT, file_id INT, NT_sequence LONGTEXT, AA_sequence TEXT,  description TEXT, type char(2), length INT, date DATE)ENGINE=innodb";
  $dbh->do($tableStr);
}

#create the directory for the database file, if it doesn't exist
if(!(-d "$dbPath\\$tablename")){
  mkdir("$dbPath\\$tablename");
}

#get the largest sb number thus far
my $sbStr = "SELECT max(end) as largest FROM sb_catalog";
$sth = $dbh->prepare($sbStr);
$sth->execute();
my $curSBID;
my $sbStart;
while(my $hashref = $sth->fetchrow_hashref()){
  $curSBID = $hashref->{largest} + 1;
  $sbStart = $curSBID;
}

#open database file
open(OUTFILE, ">>$dbPath\\$tablename\\$type" . "DatabaseFile.fas") or die "Coult not open NT database file";
#open temp load file for sequences

my $seqio = new Bio::SeqIO(-format=>"fasta", -file=>"$in");
#insert the sequences into the sequence table for the project, and write sequences to database file
my $numSeqs = 0;
while(my $seq = $seqio->next_seq){
  my $thisSBID = $curSBID;
  my $sequence = $seq->seq;
  #add the SBID to the description, if it doesn't have it already
  my $description;
  my $str;
  if($seq->display_id =~ /sb\|/){
      #if this sequence has an SBID, it's already been added to the DB, we're just updating either the NT or AA sequence
      #strip out the SBID
      $description = $seq->display_id . " " . $seq->desc;
      my @splitDesc = split(/\|/, $seq);
      #if the seq already has an SBID, then just update the sequence
      $str = "UPDATE " . $fullname . " SET " . $type . "_sequence='" . $sequence . "' WHERE sb_id='" . $splitDesc[1] . "'";
  }
  else{
      $description = "sb\|" . $thisSBID . "\| " . $seq->display_id . " " . $seq->desc;
      $description =~ s/(\"|\')//g;
      $str = "INSERT INTO " . $fullname . " (project_id, file_id, sb_id, $type\_sequence, description, length, type, date) VALUES ('" . $tablename . "', '" . $fileID . "', '" . $thisSBID . "', '". $sequence . "', '" . $description . "', '" . $seq->length . "', '" . $type . "', CURDATE())";
      $curSBID++;
  }
  print OUTFILE ">$description\n";
  print OUTFILE "$sequence\n";
  $dbh->do($str);
  $numSeqs++;
}
#add entry to the sb_catalog
my $query = "INSERT INTO sb_catalog (begin, end, fileID, projectID) VALUES ('" . $sbStart . "', '" . ($curSBID - 1) . "', '" . $fileID . "', '". $tablename . "')";
$dbh->do($query);

#call formatdb on the database file
if($type eq "NT"){
   system("C:/454/standaloneblast/bin/formatdb -i $dbPath\\$tablename\\$type" . "DatabaseFile.fas -p F");
}
elsif($type eq "AA"){
   system("C:/454/standaloneblast/bin/formatdb -i $dbPath\\$tablename\\$type" . "DatabaseFile.fas -p T");
}

#update the number of sequences for the project and the file
#get the number of seqs currently in the project
my $str;
my $mod;
if($type eq "NT"){
  $mod = "_NT_";
}
else{
  $mod = "_AA_";
}
$str = "SELECT num" . $mod . "seqs FROM project_directory WHERE projectID = '" . $tablename . "'";
$sth = $dbh->prepare($str);
$sth->execute();
my $currentCount;
while(my $hashref = $sth->fetchrow_hashref()){
  $currentCount = $hashref->{"num" . $mod . "seqs"};
}

my $newCount = $currentCount + $numSeqs;

if($type eq "NT"){
  $str = "UPDATE project_directory SET num_NT_seqs ='" . $newCount . "' WHERE projectID ='" . $tablename . "'";
  $dbh->do($str);
  $str = "UPDATE project_files SET num_seqs ='" . $numSeqs . "' WHERE fileID ='" . $fileID . "'";
}
else{
  $str = "UPDATE project_directory SET num_AA_seqs ='" . $newCount . "' WHERE projectID ='" . $tablename . "'";
  $dbh->do($str);
  $str = "UPDATE project_files SET num_seqs ='" . $numSeqs . "' WHERE fileID ='" . $fileID . "'";
}
$dbh->do($str);

