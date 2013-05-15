use strict;
use warnings;
use DBI;
use threads;
use Bio::SeqIO;
use Net::SMTP;

my $queryPath = "C:\\wamp\\www\\slimebase2\\results";
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
    my $sessionID = $hashref->{'session_id'};
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
    my $database = "";
    my $query = "$queryPath\\$folder\\query.fas";
    if(($program eq "blastn") || ($program eq "tblastn")){
      #BLAST against a nucleotide database
      $database = "$dbPath\\$resourceID\\nt" . "DatabaseFile.fas";
    }
    else{
      #BLAST against a protein database
      $database = "$dbPath\\$resourceID\\aa" . "DatabaseFile.fas";
    }
    print "Query: $query\n";
    print "Database: $database\n";
    my $outFile = "$queryPath\\$folder\\BLAST_results\\blastOutput.out";
    #run the BLAST
    print "$program started\n";
    system("C:/454/standaloneblast/blast-2.2.21+/bin/$program -db $database -query $query  -out $outFile -num_descriptions $numRes -num_alignments $numRes -evalue $eval");
    print "Counting hits\n";
    #count the results
    system("perl C:/PerlScripts/Statistics/countHitsDirection.pl $outFile $countFile $eval .01 ");
    print "Splitting BLAST Results\n";
    #split the results
    print "Drawing Graphics\n";
    #system("perl C:/PerlScripts/Parse/splitBlast.pl $outFile");
    #render the results
    system("perl C:/PerlScripts/Graphics/renderBlast2.pl $outFile $graphicsDir");
    #add entry to the finished_jobs table
    my $insert = "INSERT INTO finished_jobs (session_id, blast_id, description, db) VALUES ('" . $sessionID . "', '" . $folder . "', '" . $program . " against " . $projectName . "', '" . $resourceID . "')";
    #create a thread to pasrse the output
    my $thr = threads->create('parseBlast', $outFile, $resourceID, $sessionID, $folder);
    $thr->detach();
    $dbh->do($insert);
    if($email){
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
      $smtp->datasend("Your BLAST job has finished. You can download your results by following this link: http://150.176.130.196:8888/seq_view/results/$folder/$folder.zip");
      $smtp->dataend();
      $smtp->quit();
    }
    #remove the BLAST job from the queue
    $str = "DELETE FROM blast_queue WHERE blast_id ='" . $blastID . "'";
    $sth2 = $dbh->prepare($str);
    $sth2->execute();
  }
  #wake up every second
  sleep(1);
}

sub parseBlast{
  my ($blastOutput, $db, $sessionID, $blastID) = @_;
  my($directory, $filename) = $blastOutput =~ m/(.*\\)(.*)$/;
  open(INFILE, "<$blastOutput") or die "Could not open blast input file to parse\n";
  open(QUERIES, ">$directory\/queries_parsed.txt") or die "Could not open output file\n";
  open(HITS, ">$directory\/hits_parsed.txt") or die "Could not open output file\n";
  my @thisResult = ();
  #structure for hit array
  #0: query description
  #1: hit description
  #2: hit eval
  #3: hit identity
  #4: strand
  #5: hit sb number
  my @thisHit = ();
  my $resultID = 0;
  my $hitIndex = 0;
  my $state = "START";
  my $queryDesc = "";
  my $hitLongText = "";
  my $hitDescription = "";
  my $topDescription = "";
  my $topEval = -1;
  my $firstHitFound = 0;
  while(<INFILE>){
    chomp();
    my $line = $_;
    #remove quotes from the line
    $line =~ s/(\'|\")//g;
    if($state eq "START"){
      #if we get to this line, this is a new result
      if($line =~ /^Query=/){
        #reset result arrays
        @thisResult = ();
        $state = "QUERY";
        $line =~ s/^Query=\s+//;
        $queryDesc = $line;
        $resultID++;
        $hitIndex = 0;
      }
    }
    elsif($state eq "QUERY"){
      if($line =~ /^Length=/){
        #end of the query description, change the state to HIT analysis
        $topDescription = "";
        $topEval = -1;
        $firstHitFound = 0;
        $state = "HIT";
      }
      else{
        #append to the query description
        $queryDesc = $queryDesc . " " . $line;
      }
    }
    elsif($state eq "HIT"){

      if($line =~ /^>/){
        $hitIndex++;
        #reset the hit array
        @thisHit = ();
        $thisHit[0] = $queryDesc;
        #parse the first line of the first hit
        $hitLongText = $line;
        #strip off the greater than sign
        $line =~ s/^>//;
        $hitDescription = $line;
        #get the sb# for this hit
        ($thisHit[5]) = $line =~ m/(?<=\|)(\d+)(?=\|)/;
        $state = "HIT_DESCRIPTION";

      }
      elsif($line =~ /^\*{5} No hits found/){
        #no hits for this query sequence
       print QUERIES "INSERT INTO active_queries VALUES('" . $sessionID . "', '" . $blastID . "', '" . $resultID . "', '" . $queryDesc . "', '-1', 'N/A', '" .  $db . "', '" . 0 . "')\n";
       $state = "START";
      }
    }
    elsif($state eq "HIT_DESCRIPTION"){
      $hitLongText = $hitLongText . "<br/>" . $line;
      if($line =~ /\sScore\s\=/){
        #work some magic on the hit description
        my @splitDesc = split(/\|/, $hitDescription);
        my $tmp = "";
        for(my $i = 0; $i < scalar(@splitDesc); $i++){
          if($i == 1){
            $tmp = $tmp . "<a href=\"seqDetail.php?sbid=" . $splitDesc[1] . "&projectID=" . $db . "\" target=\"_blank\">" . $splitDesc[1] . "</a>\|";
          }
          else{
            $tmp = $tmp . $splitDesc[$i] . "|";
          }
        }
        $hitDescription = $tmp;
        #save the hit description
        $thisHit[1] = $hitDescription;
        #get the e-value
        my @split = split(/\=\s/, $line);
        $thisHit[2] = $split[2];
        #check if this is the first hit, if so, save as top
        if($firstHitFound == 0){
          $firstHitFound = 1;
          $topDescription = $hitDescription;
          $topEval = $thisHit[2];
        }
        $state = "GOT_SCORE";
      }
      else{
        $hitDescription = $hitDescription . " " . $line;
      }
    }
    elsif($state eq "GOT_SCORE"){

      $hitLongText = $hitLongText . "<br/>" . $line;
      #parse the identity
      ($thisHit[3]) = $line =~ m/(\(\d{2,3}\%\))/;
      #strip out the '(' '%' ')'
      $thisHit[3] =~ s/(\(|\%|\))//g;
      $state = "GOT_IDENT";
    }
    elsif($state eq "GOT_IDENT"){

      $hitLongText = $hitLongText . "<br/>" . $line;
      #parse the strand
      if($line =~ /Plus\/Minus/){
        $thisHit[4] = 'A';
      }
      else{
        $thisHit[4] = 'S';
      }
      $state = "ALIGNMENTS";
    }
    elsif($state eq "ALIGNMENTS"){

      $hitLongText = $hitLongText . "<br/>" . $line;
      if($line =~ /^>/){
        #write the previous hit's insert commant
        print HITS "INSERT INTO active_hits VALUES('" . $blastID . "', '" . $resultID . "', '" . $hitIndex . "', '" . $sessionID . "', '" . $thisHit[5] . "', '" . $thisHit[1] . "', '" . $thisHit[2] . "', '" . $hitLongText . "', '" . $thisHit[4] . "', '" . $thisHit[3] . "')\n";
        $hitIndex++;
        #reset the hit array
        @thisHit = ();
        $thisHit[0] = $queryDesc;
        #parse the first line of the first hit
        $hitLongText = $line;
        #strip off the greater than sign
        $line =~ s/^>//;
        $hitDescription = $line;
        #get the sb# for this hit
        ($thisHit[5]) = $line =~ m/(?<=\|)(\d+)(?=\|)/;
        $state = "HIT_DESCRIPTION";

      }
      elsif($line =~ /^Effective search space/){

        #write insert command for the query
        print QUERIES "INSERT INTO active_queries VALUES('" . $sessionID . "', '" . $blastID . "', '" . $resultID . "', '" . $queryDesc . "', '$topEval', '$topDescription', '" . $db . "', '" . $hitIndex . "')\n";
        #write the insert command for the last hit
        print QUERIES "INSERT INTO active_hits VALUES('" . $blastID . "', '" . $resultID . "', '" . $hitIndex . "', '" . $sessionID . "', '" . $thisHit[5] . "', '" . $thisHit[1] . "', '" . $thisHit[2] . "', '" . $hitLongText . "', '" . $thisHit[4] . "', '" . $thisHit[3] . "')\n";
        $state = "START";
      }

    }


  }
}
