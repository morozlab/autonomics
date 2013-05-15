<?php

session_start();
include('../../includes/dbConnect.php');
include('../../includes/config.php');

mysql_select_db('moroz_lab');

//make sure this user has access to this BLAST result
if(isset($_GET['blastID'])){
	$blastID = $_GET['blastID'];	
}
else{
	$blastID = 0;
}

$query = "SELECT * FROM finished_jobs WHERE session_id='" . session_id() . "' AND blast_id='" . $blastID . "'";
$result = mysql_query($query);
if(mysql_num_rows($result)){
	$row = mysql_fetch_array($result);
	if($row['hits_loaded'] == 0){
		$query = 'LOAD DATA LOCAL INFILE "' . $resultPath . '/' . $blastID . '/BLAST_results/hits_parsed.txt" INTO TABLE active_hits FIELDS ENCLOSED BY "\'" LINES TERMINATED BY "\tENDLINE"';
		mysql_query($query);
		$query = "UPDATE finished_jobs SET hits_loaded='1' WHERE blast_id='" . $blastID . "' AND session_id ='" . session_id() . "'";
		mysql_query($query);
	}
	$toEcho = $query;
	$query = "SELECT default_type FROM project_directory WHERE projectID='" . $row['db'] . "'";
	$defaultResult = mysql_query($query);
	$defaultRow = mysql_fetch_array($defaultResult);
	$defaultType = $defaultRow['default_type'];
	header('Content-type: text/plain');
	header('Content-Disposition: attachment; filename="' . $blastID . '"_blast_result.txt"'); 
	//create the appropriate join to get all the data I'll need
	$tableBuilder = "SELECT query_description, hit_sb, hit_description, hit_strand, hit_eval, hit_identity, " . $defaultType . "_sequence, abundance, num_hits FROM active_queries, active_hits, " . $row['db'] . "_sequences WHERE active_queries.blast_id = active_hits.blast_id AND active_queries.result_id=active_hits.result_id AND sb_id = hit_sb AND active_queries.blast_id='" . $blastID . "'";
	$tableResult = mysql_query($tableBuilder);
	echo("Query\tHit ID\tHit's Best Annotation\tHit Sequence\tAlignment E-Value\tAlignment % Identity\tHit Sequence Abundance\n");
	while($tableRow = mysql_fetch_array($tableResult)){
		if($tableRow['num_hits'] == 0){
			echo($tableRow['query_description'] . "\tNo sigifnicant alignments\n");
		}
		else{
			echo($tableRow['query_description'] . "\t" . $tableRow['hit_sb'] . "\t" . $tableRow['hit_description'] . "\t" . $tableRow[$defaultType . '_sequence'] . "\t" .  $tableRow['hit_eval'] . "\t" . $tableRow['hit_identity'] . "\t" . $tableRow['abundance'] . "\n");
		}
	}
}
?>