<?php
include('../includes/dbConnect.php');
mysql_select_db('moroz_lab');

//check if we have homology data for this project
$query = "SELECT sb_id FROM sorted_homology WHERE project_id ='" . $projectID . "' AND sort_id='1' LIMIT 0,1";
$result = mysql_query($query);
$haveAnnotation = false;
if(mysql_num_rows($result)){
	$haveAnnotation = true;
}

//check if we have GO annotation data for this project
$query = "SELECT sb_id FROM go_annotation WHERE project_id ='" . $projectID . "' LIMIT 0,1";
$result = mysql_query($query);
$haveGO = false;
if(mysql_num_rows($result)){
	$haveGO = true;
}

//check if there's KEGG annotation
$query = "SELECT has_kegg FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$haveKEGG = false;
$row = mysql_fetch_array($result);
if($row['has_kegg'] == 'Y'){
	$haveKEGG = true;	
}


//check if we have abundance for this project
$query = "SELECT has_abundance FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$haveAbundance = false;
$row = mysql_fetch_array($result);
if($row['has_abundance'] == 'Y'){
	$haveAbundance = true;
}
if(($haveAbundance == false) && ($sort == 'abundance')){
	$sort = 'evalue';
	$_SESSION['sort'] = 'evalue';
}

?>