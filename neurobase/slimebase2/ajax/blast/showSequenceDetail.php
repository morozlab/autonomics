<?php
include('../../../includes/dbConnect.php');

mysql_select_db('moroz_lab');

$seqID = $_GET['seqID'];
$projectID = $_GET['projectID'];
$program = $_GET['program'];
if($program == "blastn" || $program == "tblastn" || $program == "tblastx"){
	$mod = "NT_sequence";
}
else{
	$mod = "AA_sequence";
}
//get the sequence and the original descriptor
$query = "SELECT " . $mod . ", description FROM " . $projectID . "_sequences WHERE sb_id ='" . $seqID . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$description = ">" . $row['description'];
	$sequence = $row[$mod];
}

$sequence = chunk_split($sequence, 100);

$query = "SELECT text FROM sorted_homology, annotation_db WHERE annotation_db.annotation_id = sorted_homology.annotation_id AND sorted_homology.sb_id ='" . $seqID . "' AND sorted_homology.project_id='" . $projectID . "' AND sorted_homology.sort_id='1' ORDER BY evalue LIMIT 0,1";
$result = mysql_query($query);
$hasAnnotation = false;
while($row = mysql_fetch_array($result)){
	$hasAnnotation = true;
	$description = $description . $row['text'];
}
if($hasAnnotation == false){
	$description = $description . " No Annotation";
}
$toEcho = "<a href='seqDetail.php?sbid=" . $seqID . "&projectID=" . $projectID . "' title='view full sequence detail'>";
$toEcho .= $description ;
$toEcho .= "</a><br/>" . $sequence;
echo($toEcho);

?>