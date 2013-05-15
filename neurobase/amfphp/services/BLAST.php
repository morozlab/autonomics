<?php
include('../../includes/dbConnect.php');
include('../../includes/programPaths.php');

class BLAST{

	public function __construct(){
		mysql_select_db('moroz_lab');
	}
	
	public function prepareBLAST($params){
		$seqs = $params[0];
		$seqFile = $params[1];
		$projID = $params[2];
		$eval = $params[3];
		$email = $params[4];
		$program = $params[5];
		//make a directory to hold the results of this BLAST
		$folderName = time();
		$path = 'C:\\wamp\\www\\seq_view\\results\\' . $folderName;
		mkdir($path);
		//if there was an uploaded file, change it into the temporary query file
		if($seqFile != "none"){
			//move the temp file to the result directory
			rename('C:\\wamp\\www\\seq_view\\tmp\\' . $seqFile, $path . '\\query.fas');
		}
		else{
			//open query file for writing
			$fh = fopen($path . "\\query.fas", "ab");
			//check if the sequences came from paste box, or if we need to pull them out of DB
			if(sizeof($seqs[0]) == 2){
				//put all of the sequences in a temporary file used for the BLAST
				for($i = 0; $i < sizeof($seqs); $i++){
					fwrite($fh, ">" . $seqs[$i][0] . "\n");
					fwrite($fh, $seqs[$i][1] . "\n");
				}
			}
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
			}
			fclose($fh);
		}
		//insert the BLAST job into the BLAST queue
		$query = "INSERT INTO blast_queue SET resource_id ='" . $projID . "', folder ='" . $folderName . "', eval='" . $eval . "', email='" . $email . "', program ='" . $program . "', resource_type='project'";
		$result = mysql_query($query);
		$array = array();
		$array[0] = $folderName;
		$array[1] = mysql_insert_id();
		return $array;
		
	}
	
	function loadBlastResult($blastID){
		$basePath = "C:\\wamp\\www\\seq_view\\results";
		$secondHalf = "\\" . $blastID . "\\BLAST_results\\blastOutput.out";
		$fh = fopen($basePath . $secondHalf, 'r');
		$returnText = array();
		while($line = fgets($fh)){
			array_push($returnText, chop($line));
		}
		fclose($fh);
		return $returnText; 
	}
	
	function loadQuantification($blastID){
		$basePath = "C:\\wamp\\www\\seq_view\\results";
		$secondHalf = "\\" . $blastID . "\\Quantification\\blastOutputQuantification.txt";
		$fh = fopen($basePath . $secondHalf, 'r');
		$returnText = array();
		while($line = fgets($fh)){
			array_push($returnText, chop($line));
		}
		return $returnText; 
	}
	
	function getGraphicFilenames($blastID){
		$basePath = "C:\\wamp\\www\\seq_view\\results";
		$secondHalf = "\\" . $blastID . "\\Graphics";
		$dh = opendir($basePath . $secondHalf);
		$returnArray = array();
		while($filename = readdir($dh)){
			if(preg_match('/(^\.$|^\.\.$)/', $filename)){
				continue;
			}
			array_push($returnArray, $filename);
		}
		return $returnArray;
	}
	
	function downloadCustom($seqArray){
		//create a temporary filename
		$fileName = "custom" . time() . ".fas";
		//while the file exists, change the name
		while(file_exists("../../seq_view/tmp/" . $fileName)){
			$fileName = "custom" . $time() . "\.fas";
		}
		$fh = fopen("../../seq_view/tmp/" . $fileName, 'w');
		//get the sequence for each sequence in the array
		for($i = 0; $i < sizeof($seqArray); $i++){
			$proj = $seqArray[$i][1];
			$seqID = $seqArray[$i][0];
			$type = $seqArray[$i][2];
			$query = "SELECT " . $type . "_sequence, description FROM " . $proj . "_sequences WHERE seq_id ='" . $seqID . "'";
			$result = mysql_query($query);
			while($row = mysql_fetch_array($result)){
				if(is_null($row[$type . "_sequence"])){
					//need to translate the sequence
					//first, get the NT sequence
					$query2 = "SELECT NT_sequence FROM " . $proj . "_sequences WHERE seq_id ='" . $seqID . "'";
					$result2 = mysql_query($query2);
					while($row2 = mysql_fetch_array($result2)){
						$ntSeq = $row2['NT_sequence'];
					}
					$array = array();
					exec("perl C:\PerlScripts\Format\webTranslate.pl " . $ntSeq, $array);
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
					fwrite($fh, $row[$type . "_sequence"] . "\n");
				}
			}	
		}
		fclose($fh);
		$command = "C:\\7-Zip\\7z" . " a -r -tzip " . "C:\\wamp\\www\\seq_view\\tmp\\" . $fileName . ".zip " . "C:\\wamp\\www\\seq_view\\tmp\\" . $fileName;
		$return = array();
		exec($command, $return);
		//remove the original file
		unlink("../../seq_view/tmp/" . $fileName);
		//return the name of the zipped file
		return($fileName . ".zip");
		
	}
	
	public function doNothing(){
	
	}

}


?>