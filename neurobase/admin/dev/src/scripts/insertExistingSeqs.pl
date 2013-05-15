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
my $dbPath = "C:\\wamp\\www\\seq_view\\database";
my $fullname = $tablename . "_sequences";

#open database file
open(OUTFILE, ">>$dbPath\\$tablename\\$type" . "DatabaseFile.fas") or die "Coult not open database file";
my $seqio = new Bio::SeqIO(-format=>"fasta", -file=>"$in");
#insert the sequences into the sequence table for the project, and write sequences to database file
my $numSeqs = 0;
while(my $seq = $seqio->next_seq){
  my $sequence = $seq->seq;
  #add the SBID to the description, if it doesn't have it already
  my $description = $seq->display_id . " " . $seq->desc;
  my $str = "UPDATE " . $fullname . " SET " . $type . "_sequence='" . $sequence . "' WHERE description LIKE 'sb\|%\| " . $description . "'";
  print OUTFILE ">$description\n";
  print OUTFILE "$sequence\n";
  $dbh->do($str);
  $numSeqs++;
}

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
my $sth = $dbh->prepare($str);
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

