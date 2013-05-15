<?php

include('../../includes/dbConnect.php');
$projectID = $_GET['projectID'];
if(!(is_numeric($projectID))){
	die;
}
//select database
mysql_select_db('moroz_lab');


//for the specified project, determine the number of P, F, and C GO entries and return them as an array
$resultString = "";
$array = array(6490, 2668, 4396);
foreach($array as $compartment){
	/*$query = "SELECT COUNT(go_annotation.sb_id) as c FROM go_annotation, go_catalog WHERE ((go_annotation.go_higher = go_catalog.go_id) OR (go_annotation.go_higher = 'NA' and go_catalog.go_id = go_annotation.go_id) OR (go_annotation.go_id = go_catalog.go_id and go_annotation.go_higher NOT IN (SELECT go_id FROM go_catalog WHERE compartment='" . $compartment . "'))) AND go_catalog.compartment ='" . $compartment . "' AND go_annotation.project_id ='" . $projectID . "'";
	
	$query = "SELECT COUNT(go_annotation.sb_id) as c FROM go_annotation, go_catalog WHERE go_annotation.go_higher = go_catalog.go_id and go_annotation.project_id ='" . $projectID . "' AND go_catalog.compartment ='" . $compartment . "'";
	$count += addToCount($query);
	$query = "SELECT COUNT(go_annotation.sb_id) as c FROM go_annotation, go_catalog WHERE go_annotation.go_higher ='NA' and go_catalog.go_id = go_annotation.go_id and go_annotation.project_id = '" . $projectID . "' and go_catalog.compartment ='" . $compartment . "'";
	$count += addToCount($query);
	$query = "SELECT COUNT(go_annotation.sb_id) as c FROM go_annotation, go_catalog WHERE (go_annotation.go_higher NOT IN (SELECT DISTINCT(go_id) FROM go_catalog WHERE compartment='" . $compartment . "')) AND go_annotation.go_higher !='NA' AND go_catalog.go_id = go_annotation.go_id AND go_annotation.project_id ='" . $projectID . "' AND go_catalog.compartment ='" . $compartment ."'";
	$count += addToCount($query);*/
	$query = "SELECT * FROM go_categories WHERE project_id='" . $projectID . "' AND id='" . $compartment . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$resultString = $resultString . "<div style='padding-bottom: 4px; width:880px; float:left;' id='" . $compartment . "_div'" . '><a id="' . $compartment . '" target="_self" href="#" onclick="toggleCategory(this.id, \'expand\', \'top_level\');return(false);">' . $row['name'] . " (" . $row['unique_seqs'] . " | " . $row['unique_annotations'] . ')</a><div id="' . $compartment . "_children" . '"></div></div>';
}

echo($resultString);

function addToCount($query){
	$result = mysql_query($query);
	$count = 0;
	while($row = mysql_fetch_array($result)){
		//for each top-level, build the div structure needed to display things
		$count += $row['c'];
	}
	return $count;
}


?>