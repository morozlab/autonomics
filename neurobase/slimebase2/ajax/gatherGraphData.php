<?php

include('../../includes/dbConnect.php');

$projectID = $_GET['projectID'];
$level = $_GET['level'];

mysql_select_db('moroz_lab');
$check = "SELECT * FROM go_categories WHERE project_id='" . $projectID . "' LIMIT 0,1";
$checkResult = mysql_query($check);
$resultString = "";
if(!(mysql_num_rows($checkResult))){
if($level == 'top'){
	

	#count the number of each of the major categories for the project
	$array = array('P', 'F', 'C');

	$resultArray = array();
	$index = 0;
	foreach($array as $compartment){
		$query = "SELECT COUNT(go_annotation.sb_id) as num FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_annotation.go_id = go_catalog.go_id and go_catalog.compartment ='" . $compartment . "'";
		$result = mysql_query($query);
		while($row = mysql_fetch_array($result)){
			$resultArray[$index] = $row['num'];
		}
		$index++;
	}
	echo(json_encode($resultArray));
}

else{
	$returnArray = array();
	$returnArray = countCategories($level);
	echo(json_encode($returnArray));
}
}
else{
	$array = array(6490, 2668, 4396);
	$resultArray = array();
	$index = 0;
	foreach($array as $term){
		$query = "SELECT unique_annotations FROM go_categories WHERE project_id='" . $projectID . "' AND id='" . $term . "'";
		$result = mysql_query($query);
		while($row = mysql_fetch_array($result)){
			$resultArray[$index] = $row['unique_annotations'];	
		}
		$index++;
	}
	echo(json_encode($resultArray));
}


function countCategories($level){
	global $projectID;
	$returnArray = array();
	$query = "SELECT DISTINCT(go_catalog.go_id), go_catalog.description FROM go_catalog, go_annotation WHERE go_catalog.go_id = go_annotation.go_higher AND go_annotation.project_id ='" . $projectID . "' AND go_catalog.compartment='" . $level . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		//count the sequences for each id
		$query = "SELECT COUNT(sb_id) as num FROM go_annotation WHERE project_id ='" . $projectID . "' AND go_higher ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		while($row2 = mysql_fetch_array($result2)){
			$count = $row2['num'];
		}
		$returnArray[$row['description']] = $count;
	}
	$query = "SELECT DISTINCT(go_catalog.go_id), go_catalog.description FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_catalog.go_id = go_annotation.go_id AND go_catalog.compartment ='" . $level . "' AND go_annotation.go_higher='NA' ORDER BY go_catalog.description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as num FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$slimCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$count = $row2['num'];
		}
		$returnArray[$row['description']] = $count;
	}
	
	/*$query = "SELECT DISTINCT(go_catalog.go_id), go_catalog.description FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_catalog.go_id = go_annotation.go_id AND go_catalog.compartment ='" . $level . "' AND go_annotation.go_higher !='NA' AND go_annotation.go_higher NOT IN (SELECT DISTINCT(go_id) FROM go_catalog WHERE compartment='" . $level . "') ORDER BY go_catalog.description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as num FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$slimCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$count = $row2['num'];
		}
		$returnArray[$row['description']] = $count;
	}*/
	return $returnArray;
}
?>