<?php
include('../../includes/dbConnect.php');

$projectID = $_GET['projectID'];
$extra = $_GET['extra'];
$term = $_GET['term'];

mysql_select_db('moroz_lab');


$returnString = "";

//first, get all the pertinent information for the children of this term
$query = "SELECT term.acc, term.name, unique_seqs, unique_annotations, term.id FROM term, term2term, go_categories WHERE term2term.term1_id='" . $term . "' AND term.id = term2_id AND go_categories.project_id='" . $projectID . "' AND go_categories.id = term2term.term2_id ORDER BY term.name";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$goID = $row['id'];
	$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $goID . "_div'><a href='#' target='_self' id='" . $goID . "' onclick='toggleCategory(this.id, \"expand\", \"slim\");return(false);'>" . $row['name'] . " (" . $row['unique_seqs'] . " | " . $row['unique_annotations'] . ")</a><div id='" . $goID . "_children'></div></div>"; 	
}

//get all of the sequences falling into this category
$query = "SELECT acc FROM term WHERE id='" . $term . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$GO = $row['acc'];
$array = preg_split("/\:/", $GO);
$GO = $array[1];
$query = "SELECT go_annotation.sb_id, go_annotation.go_id FROM go_annotation WHERE go_annotation.project_id = '" . $projectID . "' AND go_annotation.go_id ='" . $GO . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$sb = $row['sb_id'];
	$query2 = "SELECT text, evalue FROM homology, annotation_db WHERE homology.sb_id = '" . $sb ."' AND homology.annotation_id = annotation_db.annotation_id ORDER BY evalue LIMIT 0,1";
	$result2 = mysql_query($query2);
	$found = false;
	while($row2 = mysql_fetch_array($result2)){
		$found = true;
	$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| " . $row2['text'] . " evalue: " . $row2['evalue'] . "</a></div>"; 
	}
	if($found == false){
		$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| No Annotation</a></div>";
	}
}

echo($returnString);

/*
if($extra == 'top_level'){
	if($target == "P"){
		$id = 6490;	
	}
	else if($target == "C"){
		$id = 4396;
	}
	else if($target == "F"){	
		$id = 2668;
	}
	
	$query = "SELECT acc, name FROM term2term, term WHERE term.id = term2_id AND term1_id='" . $id . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$array = preg_split("/\:/", $row['acc']);
		if(sizeof($array) != 1){
			$goID = $array[1];
			$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $goID . "_div'><a href='#' target='_self' id='" . $goID . "' onclick='toggleView(this.id, \"expand\", \"slim\");return(false);'>" . $row['name'] . "</a><div id='" . $goID . "_children'></div></div>"; 
		}
	}
	
	
}
else if($level == 'slim'){
	$query = "SELECT DISTINCT(go_catalog.go_id), description FROM go_catalog, go_annotation WHERE go_annotation.go_higher = '" . $term ."' AND go_annotation.go_id = go_catalog.go_id and go_annotation.project_id ='" . $projectID . "' ORDER BY description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as num FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$seqCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$seqCount = $row2['num'];
		}
		$returnString .= "<div style='float:left; width:840px; padding-bottom:2px; padding-left: 40px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"sequence\"); return(false);'>" . $row['description'] . " (" . $seqCount . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 
	}
}

else if($level == 'sequence'){
	$query = "SELECT go_annotation.sb_id, go_annotation.go_id FROM go_annotation WHERE go_annotation.project_id = '" . $projectID . "' AND go_annotation.go_id ='" . $term . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$sb = $row['sb_id'];
		$query2 = "SELECT text, evalue FROM homology, annotation_db WHERE homology.sb_id = '" . $sb ."' AND homology.annotation_id = annotation_db.annotation_id ORDER BY evalue LIMIT 0,1";
		$result2 = mysql_query($query2);
		$found = false;
		while($row2 = mysql_fetch_array($result2)){
			$found = true;
		$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| " . $row2['text'] . " evalue: " . $row2['evalue'] . "</a></div>"; 
		}
		if($found == false){
			$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| No Annotation</a></div>";
		}
	}
}
echo($returnString);*/


?>