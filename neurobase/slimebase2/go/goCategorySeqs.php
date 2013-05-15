<?php
session_start();
include('../../includes/dbConnect.php');
include('../../includes/restrictAccess.php');

$projectID = $_GET['projectID'];
if(isset($_GET['level'])){
	$level = $_GET['level'];
}
$term = $_GET['term'];
if(isset($_GET['expected'])){
	$expected = $_GET['expected'];
}

$returnString = "";
//get all of the sequences falling into this category
$query = "SELECT acc FROM term WHERE id='" . $term . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$GO = $row['acc'];
$array = preg_split("/\:/", $GO);
$GO = $array[1];
$query = "SELECT go_annotation_new.sb_id, go_annotation_new.go_id FROM go_annotation_new, " . $projectID . "_sequences WHERE go_annotation_new.sb_id = " . $projectID . "_sequences.sb_id AND go_annotation_new.go_id ='" . $GO . "'";
	$result = mysql_query($query);	

while($row = mysql_fetch_array($result)){
	$sb = $row['sb_id'];
	$query2 = "SELECT annot, eval FROM best_annotations WHERE sb_id = " . $sb ." AND project_id=" . $projectID;
	$result2 = mysql_query($query2);
	$found = false;
	while($row2 = mysql_fetch_array($result2)){
		$found = true;
		$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a style='color:#063480;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| " . $row2['annot'] . " evalue: " . $row2['eval'] . "</a></div>"; 
	}
	if($found == false){
		$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#063480;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| No Annotation</a></div>";
	}
}
echo($returnString);

?>