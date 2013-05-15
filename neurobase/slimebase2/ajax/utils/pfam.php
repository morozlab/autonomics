<?php

session_start();
include('/var/www/html/neurobase/includes/restrictAccess.php');
include('/var/www/html/neurobase/includes/dbConnect.php');

$method = $_GET['method'];
$args = $_GET['args'];
$argArray = explode(",", $args);
//first argument is always the projectID
$projectID = $argArray[0];
//check permissions
include('/var/www/html/neurobase/includes/slimebase2/security/checkPermissions.php');
	

if($method == "getPfamDetails"){
	$pfamID = $argArray[1];
	$query = "SELECT pfamA_acc, description, comment FROM pfama WHERE pfamA_acc='" . $pfamID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	//add the pfam description and close button
	$toEcho = "<div id='pfam-detail-header'>" . $row['description'] . "(" . $row['pfamA_acc'] . ")</div><div id='pfamDetailClose'><a href='#' onclick='closePfamDetail();return(false);'>close x</a></div>";
	//add the full HTML info from Pfam (later)
	//add the comment from this entry
	$toEcho .= "<div id='pfam-detail-comment'>" . $row['comment'] . "</div>";
	//get the sequences from this project with this pfam domain
	$toEcho .= "<div id='pfam-sequences-header'>Sequences with '" . $row['description'] . "' domain: </div>";
	$query = "SELECT pfam_annotations.sb_id, description, length, evalue FROM " . $projectID . "_sequences, pfam_annotations WHERE " . $projectID . "_sequences.sb_id = pfam_annotations.sb_id AND pfam_annotations.project_id='" . $projectID . "' AND pfam_annotations.pfamA_acc='" . $pfamID . "' ORDER BY evalue";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		//draw the pfam domains for each sequence
		$toEcho .= "<div class='pfam-detail-diagram-heading'>" . $row['sb_id'] . " " . $row['description'] . " (" . $row['evalue'] . ")</div>";
		$toEcho .= "<div id='pfam-detail-" . $row['sb_id'] . "-rendering'></div>";
		$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, description, evalue, start, stop, type FROM pfam_annotations, pfama WHERE project_id='" . $projectID . "' AND sb_id='" . $row['sb_id'] . "' AND pfam_annotations.pfamA_acc = pfama.pfamA_acc";	
		$pfamDomainResult = mysql_query($query);
		$toEcho .= "<script language='javascript' type='text/javascript'>featureArray = new Array();";
		while($row2 = mysql_fetch_array($pfamDomainResult)){
			$toEcho .= "featureArray.push(new pfamAnnotation('" . $row2['pfamA_acc'] . "', '" . $row2['start'] . "', '" . $row2['stop'] . "', '" . $row2['description'] . "', '" . $row2['type'] . "', '" . $row2['pfamA_id'] . "'));\n";	
		}
		$toEcho .= "displayPfamGraphic(featureArray, '" . $row['length'] . "', '" . $row['sb_id'] . "');";
		$toEcho .= "</script>";
	}
	echo($toEcho);
}
else if($method == 'loadPage'){
	$argArray = explode(",", $args);
	$projectID = $argArray[0];
	$character = $argArray[1];
	$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, pfama.description, pfama.type, counts FROM pfama, pfam_domain_counts WHERE pfama.pfamA_acc = pfam_domain_counts.acc AND project_id='" . $projectID . "' AND pfama.pfamA_id LIKE '" . $character . "%'";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$toEcho = generatePfamTable($result, $projectID);
	}
	echo($toEcho);	
}
else if($method == "updateSearchRes"){
	$argArray = explode(",", $args);
	$projectID = $argArray[0];
	$search = $argArray[1];
	$search = preg_replace("/(^\s+|\s+$)/", "", $search);
	$searchArray = explode(" ", $search);
	$newSearch = "";
	if(sizeof($searchArray) == 1){
		$newSearch = "+\"" . $searchArray[0] . "\""; //. "*";	
	}
	else{
		for($i = 0; $i < sizeof($searchArray); $i++){
			if(preg_match('/\-/', $searchArray[$i])){
				$searchArray[$i] = "\"" . $searchArray[$i] . "\"";
			}
			if($i == 0){
				$newSearch .= "+" . $searchArray[$i]; 
			}
			else{
				$newSearch .= " +" . $searchArray[$i];
			}
		}
	}
	$search = $newSearch;
	$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, pfama.description, pfama.type, counts FROM pfama, pfam_domain_counts WHERE pfama.pfamA_acc = pfam_domain_counts.acc AND project_id='" . $projectID . "' AND MATCH(pfama.description) AGAINST ('" . $search . "' IN BOOLEAN MODE) LIMIT 0,200";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$toEcho = generatePfamTable($result);
	}
	else{
		$toEcho = 'null';	
	}
	echo($toEcho);
}

function generatePfamTable($result){
	global $projectID;
	$toEcho = "<table id='pfam-browsing-table'>";
	$toEcho .= "<tr class='pfam-table-header'><td>id</td><td>accession</td><td>description</td><td class='type-column'>type</td><td># of seqs</td></tr>";
	$colorator = 0;
	while($row = mysql_fetch_array($result)){
		if($colorator == 0){
			$toEcho .= "<tr class='pfam-table-light' id=\"" . $row['pfamA_acc'] . "\"><td>" . $row['pfamA_id'] . "</td><td>" . $row['pfamA_acc'] . "</td><td>" . $row['description'] . "</td><td>" . $row['type'] . "</td><td>" . $row['counts'] . "</td></tr>";
			$colorator = 1;
		}
		else{
			$toEcho .= "<tr class='pfam-table-dark' id=\"" . $row['pfamA_acc'] . "\"><td>" . $row['pfamA_id'] . "</td><td>" . $row['pfamA_acc'] . "</td><td>" . $row['description'] . "</td><td>" . $row['type'] . "</td><td>" . $row['counts'] . "</td></tr>";
			$colorator = 0;
		}
	}
	$toEcho .= "</table>";
	return $toEcho;
}

?>