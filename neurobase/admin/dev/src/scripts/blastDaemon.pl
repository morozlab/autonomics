use strict;
use warnings;
use DBI;
use threads;
use Bio::SeqIO;
use Net::SMTP;

my $queryPath = "C:\\wamp\\www\\seq_view\\results";
my $dbPath = "C:\\wamp\\www\\seq_view\\database";
my $numRes = '1000000';
my $dsn = "dbi:mysql:database=moroz_lab;host=localhost";
my $dbh = DBI->connect($dsn, "root", "meow12");

while(1){
  my $str = 'SELECT * FROM blast_queue LIMIT 0,1';
  my $sth = $dbh->prepare($str);
  $sth->execute();
  while(my $hashref = $sth->fetchrow_hashref()){
    #there is at least one BLAST waiting in the queue, do the BLAST
    my $blastID = $hashref->{'blast_id'};
    my $resourceID = $hashref->{'resource_id'};
    my $folder = $hashref->{'folder'};
    my $resourceType = $hashref->{'resource_type'};
    my $eval = $hashref->{'eval'};
    my $email = $hashref->{'email'};
    my $program = $hashref->{'program'};
    #get the project name of the project being BLASTed against
    my $str2 = "SELECT project_name FROM project_directory WHERE projectID ='" . $resourceID . "'";
    my $sth2 = $dbh->prepare($str2);
    $sth2->execute();
    my $projectName = "";
    while(my $resref = $sth2->fetchrow_hashref()){
      $projectName = $resref->{'project_name'};
      print "Beginning BLAST $blastID against: $projectName" . "\n";
    }
    #make the BLAST result directory
    if(!(-d "$queryPath\\$folder\\BLAST_results")){
      mkdir("$queryPath\\$folder\\BLAST_results");
    }
    mkdir("$queryPath\\$folder\\Quantification");
    mkdir("$queryPath\\$folder\\Graphics");
    my $countFile = "$queryPath\\$folder\\Quantification\\blastOutputQuantification.txt";
    my $graphicsDir = "$queryPath\\$folder\\Graphics\\";
    if(($program eq "blastn") || ($program eq "tblastn")){
      #run the BLAST of the query sequences against the database
      #determine the alphabet of the sequences
      my $seqio = new Bio::SeqIO(-format=>"fasta", -file=>"$queryPath\\$folder\\query.fas");
      my $seq = $seqio->next_seq;
      #need to run tblastn, protein vs nucleotide
      my $database = "$dbPath\\$resourceID\\nt" . "DatabaseFile.fas";
      my $query = "$queryPath\\$folder\\query.fas";
      print "Query: $query\n";
      print "Database: $database\n";
      my $outFile = "$queryPath\\$folder\\BLAST_results\\blastOutput.out";

      #run the BLAST
      print "$program started\n";
      system("C:/454/standaloneblast/bin/blastall -p $program -d $database -i $query  -o $outFile -v $numRes -b $numRes -e $eval");

      print "Counting hits\n";
      #count the results
      system("perl C:/PerlScripts/Statistics/countHitsDirection.pl $outFile $countFile $eval .01 ");
      print "Splitting BLAST Results\n";
      #split the results
      print "Drawing Graphics\n";
      #system("perl C:/PerlScripts/Parse/splitBlast.pl $outFile");
      #render the results
      system("perl C:/PerlScripts/Graphics/renderBlast2.pl $outFile $graphicsDir");
    }
    else{
      #BLAST against a protein database
      my $seqio = new Bio::SeqIO(-format=>"fasta", -file=>"$queryPath\\$folder\\query.fas");
      my $seq = $seqio->next_seq;
      #need to run tblastn, protein vs nucleotide
      my $database = "$dbPath\\$resourceID\\aa" . "DatabaseFile.fas";
      my $query = "$queryPath\\$folder\\query.fas";
      print "Query: $query\n";
      print "Database: $database\n";

      my $outFile = "$queryPath\\$folder\\BLAST_results\\blastOutput.out";
      #start the BLAST

      print "$program started\n";
      system("C:/454/standaloneblast/bin/blastall -p $program -d $database -i $query  -o $outFile -v $numRes -b $numRes -e $eval");

      print "Counting hits\n";
      #count the results
      system("perl C:/PerlScripts/Statistics/countHitsDirection.pl $outFile $countFile $eval .01");
      print "Splitting BLAST Results\n";
      #split the results
      print "Drawing Graphics\n";
      #system("perl C:/PerlScripts/Parse/splitBlast.pl $outFile");
      #render the results
      system("perl C:/PerlScripts/Graphics/renderBlast2.pl $outFile $graphicsDir");
    }
    #zip up the results
    system("C:/7-Zip/7z a -r -tzip $queryPath\\$folder\\$folder.zip $queryPath\\$folder");
    #send an email about the results
    my $smtp = Net::SMTP->new('whitney.ufl.edu');
    $smtp->auth('morozgenomics', 'mg123');
    $smtp->mail("morozgenomics\@whitney.ufl.edu");
    $smtp->to($email);
    $smtp->data();
    $smtp->datasend("Subject: BLAST Job $folder Finished\n");
    $smtp->datasend("To: $email\n");
    $smtp->datasend("\n");
    $smtp->datasend("Your BLAST job has finished. You can download your results by following this link: http://10.41.128.72/seq_view/results/$folder/$folder.zip");
    $smtp->dataend();
    $smtp->quit();
    #remove the BLAST job from the queue
    $str = "DELETE FROM blast_queue WHERE blast_id ='" . $blastID . "'";
    $sth2 = $dbh->prepare($str);
    $sth2->execute();
  }
  #wake up every second
  sleep(1);
}