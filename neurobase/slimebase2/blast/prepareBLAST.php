<?php
session_start();
include_once('../../includes/dbConnect.php');
include_once('../../includes/restrictAccess.php');

isset($_POST['queryEntry']) ? $seqs = "yes" : $seqs = "none";
$projID = $_POST['database'];
$eval = $_POST['eval'];
$email = $_POST['email'];
$program = $_POST['program'];
$mode = $_POST['mode'];
//make a directory to hold the results of this BLAST
$folderName = time();
$path = '/var/www/html/neurobase/slimebase2/results/' . $folderName;
mkdir($path);
system("chmod a+rwx -R $path");
if($_FILES['queryFile']['size']){
	$seqs = "none";
	//move the temporary file to the result folder
	$path = $path . "/query.fas";
	if(move_uploaded_file($_FILES['queryFile']['tmp_name'], $path)){
		//don't need to do anything, file moved correctly
	}
	else{
		echo("<script language='javascript'>alert('There was an error uploading your file. Please try again.');location.href='../blast.php?view=default';</script>");
	}		
}
else if($seqs != "none"){
	//open query file for writing
	$fh = fopen($path . "/query.fas", "ab");
	fwrite($fh, $_POST['queryEntry']);
	/*
	//custom fasta
	else if(sizeof($seqs[0]) == 3){
		//need to look up the sequences
		for($i = 0; $i < sizeof($seqs); $i++){
			$seqID = $seqs[$i][0];
			$proj = $seqs[$i][1];
			$type = $seqs[$i][2];
			if($type == "AA"){
				$other = "NT";
			}
			else{
				$other = "AA";
			}
			$query = "SELECT " . $type . "_sequence, " . $other . "_sequence, description FROM " . $proj . "_sequences WHERE seq_id ='" . $seqID . "'";
			$result = mysql_query($query);
			while($row = mysql_fetch_array($result)){
				if(is_null($row[$type . "_sequence"])){
					if($type == "AA"){
						//need to translate the sequence
						$array = array();
						exec("perl C:\PerlScripts\Format\webTranslate.pl " . $row[$other . "_sequence"], $array);
						$seenCounter = 1;
						$identifier = $row['description'];
						$mod = "";
						for($j = 0; $j < (sizeof($array) - 1); $j++){
							if(preg_match("/^>/", $array[$j])){
								if($seenCounter == 1){
									$mod = "frame 1";
								}
								else if($seenCounter == 2){
									$mod = "frame 2";
								}
								else if($seenCounter == 3){
									$mod = "frame 3";
								}
								else if($seenCounter == 4){
									$mod = "frame -1";
								}
								else if($seenCounter == 5){
									$mod = "frame -2";
								}
								else if($seenCounter == 6){
									$mod = "frame -3";
								}
								fwrite($fh, ">" . $identifier . " " . $mod . "\n");
								$seenCounter++;
							}
							else{
								fwrite($fh, $array[$j] . "\n");
							}
						}
					}
					else{
						fwrite($fh, ">" . $row['description'] . "\n");
						fwrite($fh, $row[$other . "_sequence"] . "\n\n");
					}
				}
				else{
					fwrite($fh, ">" . $row['description'] . "\n");
					fwrite($fh, $row[$type . "_sequence"] . "\n\n");
				}
			}
		}
	}*/
	fclose($fh);
}
//insert the BLAST job into the BLAST queue
$query = "INSERT INTO blast_queue SET resource_id ='" . $projID . "', folder ='" . $folderName . "', eval='" . $eval . "', email='" . $email . "', program ='" . $program . "', resource_type='project', session_id='" . session_id() . "', mode='" . $mode . "'";
$result = mysql_query($query);
$_SESSION['view'] = 'queue';
echo("<script language='javascript' type='text/javascript'>location.href='../blast.php'</script>");
?>