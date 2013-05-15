<?php
session_start();
include('../../../includes/dbConnect.php');
include('../../../includes/config.php');

//need to perform security check and update session info here

mysql_select_db('moroz_lab');

$blastID = $_GET['blastID'];
$resultID = $_GET['resultID'];
$hitID = $_GET['hitID'];

//check if hits are loaded for this result yet
$query = "SELECT blast_id, db, hits_loaded FROM finished_jobs WHERE blast_id ='" . $blastID . "' AND session_id ='" . session_id() . "'";
$result = mysql_query($query);
//if they aren't loaded, load them
while($row = mysql_fetch_array($result)){
	$projectID = $row['db'];
	if($row['hits_loaded'] == 0){
		$query = 'LOAD DATA LOCAL INFILE "' . $resultPath . '/' . $blastID . '/BLAST_results/hits_parsed.txt" INTO TABLE active_hits FIELDS ENCLOSED BY "\'" LINES TERMINATED BY "\tENDLINE"';
		mysql_query($query) or die(mysql_error());
		$query = "UPDATE finished_jobs SET hits_loaded='1' WHERE blast_id='" . $blastID . "' AND session_id ='" . session_id() . "'";
		mysql_query($query);
	}
}

//next, get the full description for the hit and return it
$query = "SELECT hit_alignment, hit_sb FROM active_hits WHERE blast_id='" . $blastID . "' AND session_id ='" . session_id() . "' AND result_id ='" . $resultID . "' AND hit_id ='" . $hitID . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$sb = $row['hit_sb'];
	$alignment = $row['hit_alignment'];
}
$alignment .= "<br/>";
//get the hit sequence
$query = "SELECT NT_sequence, description, text, evalue FROM " . $projectID . "_sequences as t1 LEFT JOIN homology ON t1.sb_id = homology.sb_id LEFT JOIN annotation_db ON homology.annotation_id = annotation_db.annotation_id WHERE t1.sb_id ='" . $sb . "' ORDER BY evalue LIMIT 0,1";
$result = mysql_query($query);
$row = mysql_fetch_array($result); 
$exploded = preg_split("/\s+/", $row['description']);
$seqID = $exploded[0] . " ";
if(isset($exploded[1])){
	$seqID .= $exploded[1];
}
$hitDef = $seqID;
if(!(is_null($row['evalue']))){
	$seqID .=  " " . $row['text'];
	$hitDef = $seqID;
	$seqID .= " (" . $row['evalue'] . ")";
}

echo("<b>Alignment of your sequence to: </b><br/><br/><a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "' title='View Full Sequence Detail' target='_blank'>>" . $seqID . "</a><br/>");
echo(chunk_split($row['NT_sequence'], 100));
//echo(">" . $hitDef . "<br/>");
echo("<br/>" . $alignment . "<br/>");
?>