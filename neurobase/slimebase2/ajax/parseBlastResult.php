<?php

include_once('../../includes/dbConnect.php');

session_start();

mysql_select_db('moroz_lab');

$blastID = $_GET['blastID'];
$_SESSION['loadedResult'] = $blastID;
$path = "/var/www/html/neurobase/slimebase2/";
$result_path = '/slimebase2/results/' . $blastID . "/BLAST_results/blastOutput.out";

//check if we need to load the result into the active result table
$query = "SELECT session_id FROM active_queries WHERE session_id ='" . session_id() . "' AND blast_id ='" . $blastID . "'";
$result = mysql_query($query);
if(mysql_num_rows($result)){
	//there's a result that's currently loaded, just parse and display it
	parseBlastResult();
	
}
else{
	//need to load the result
	$path = $path . 'results/' . $blastID . '/BLAST_results/queries_parsed.txt';
	$query = 'LOAD DATA LOCAL INFILE "' . $path . '" INTO TABLE active_queries FIELDS ENCLOSED BY "\'" LINES TERMINATED BY "\tENDLINE"';
	mysql_query($query) or die(mysql_error());
	//result loaded, now we need to parse things from the db
	parseBlastResult();
}

function parseBlastResult(){
	global $blastID;
	global $result_path;
	$session = session_id();
	//get the program and database used in BLAST
	$query = "SELECT description FROM finished_jobs where blast_id ='" . $blastID . "' AND session_id ='" . $session . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$blastHeader = "<div id='blastDescription'><h1>" . $row['description'] . "</h1></div>";
	//get all of the results for the BLAST
	$query = "SELECT * FROM active_queries WHERE session_id ='" . $session . "' AND blast_id ='" . $blastID . "'";
	$result = mysql_query($query);
	$numResults = mysql_num_rows($result);
	$blastOptions = "<div id='blastResultOptions'><div class='blastOptionHeader'>" . $numResults . " result(s)</div>";
	$blastOptions .= "<div class='blastExportOptions'><a href='download/dlBlast.php?blastID=" . $blastID . "' target='_blank'>export result</a></div>";
	$blastOptions .= "<div class='blastExportOptions'><a href='" . $result_path . "' target='_blank'>full BLAST output</a></div";
	$queryIndex = 1;
	$toEcho = "";
	$blastOptions .= "<div class='jumpSelect'><br/><form>jump to result for: <select name='querySelect' id='querySelect' onchange='jumpResult(this.options[this.selectedIndex].value);'>"; 
	while($row = mysql_fetch_array($result)){
		//need to generate the select box at the same time we do this
		if($row['query_description'] == ""){
			$description = "Query " . $queryIndex;
		}
		else{
			$description = substr($row['query_description'], 0, 80);
		}
		$blastOptions .= "<option value='" . $queryIndex . "'>" . $description . "</option>";
		$toEcho .= "<div class='resultContainer' id='resultContainer" . $queryIndex . "'>";
		//query description header
		//start container to hold query identifier
		//$toEcho .= "<div class='innerResultBox'>";
		$toEcho .= "<div class='querySubHeader'><a name='result" . $queryIndex . "'>Query=</a></div>";
		$toEcho .= "<div class='queryDescription'>" . $description . "</div>";
		//$toEcho .= "<div class='queryHits'>" . $row['num_hits'] . "</div>";
		//$toEcho .= "</div>";
		//start container to hold the top hit info, and the other hits when the user expands
		//$toEcho .= "<div class='hitInfo'>";
		$toEcho .= "<div class='blastResultSubHeader'>Sequences producing significant alignments:</div>";
		$toEcho .= "<div class='hitSummaries'>";
		if($row['hit_summary'] == 'none'){
			$toEcho .= "No significant alignments";
		}
		else{
			$toEcho .= "<div class='summaryHeader'><div class='sbCell'>sb #</div><div class='hitDescHeader'>sequence name</div><div class='evalCell'>e-val</div><div class='idenCell'>iden%</div><div class='scoreCell'>score</div></div>" . $row['hit_summary'];
		}
		//$toEcho .= "</div>";
		//div to be target of ajax command to load hit text
		//end container for hit info
		$toEcho .= "</div>";
		//end container for result
		$toEcho .= "</div>";
		$queryIndex++;
	}
	$blastOptions .= "</select></form></div><a name='top'></div>";
	echo($blastHeader . $blastOptions . $toEcho);
}


?>