<?php
session_start();
$method = $_GET['method'];
$args = $_GET['args'];
$argArray = explode(",", $args);
//first argument is always the projectID
$projectID = $argArray[0];
$_GET['projectID'] = $projectID;
include('../../../includes/restrictAccess.php');
include('../../../includes/dbConnect.php');
mysql_select_db('moroz_lab');


//check permissions
include('../../../includes/slimebase2/security/checkPermissions.php');

if($method == "colorPath"){
	//get all of the sb_ids of the sequences matching this pathway in the project
	$pathwayID = $argArray[1];
	
	$query = "SELECT STRAIGHT_JOIN * FROM kegg_annotations JOIN best_annotations ON best_annotations.sb_id=kegg_annotations.sb_id JOIN kegg_confs ON kegg_annotations.path_id = kegg_confs.path_id WHERE kegg_annotations.project_id='" . $projectID . "' AND kegg_annotations.path_id='" . $pathwayID . "' AND kegg_annotations.project_id=best_annotations.project_id AND kegg_confs.kegg_id = kegg_annotations.kegg_id";
	$result = mysql_query($query);
	$returnArray = array();
	while($row = mysql_fetch_array($result)){
		array_push($returnArray, $row);
	}
	echo(json_encode($returnArray));
}