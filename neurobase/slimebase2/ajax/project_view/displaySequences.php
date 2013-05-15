<?php

session_start();

include_once("../../../includes/dbConnect.php");

mysql_select_db('moroz_lab');

include_once("../../../includes/slimebase2/security/ajaxCheckPermissions.php");

//check the parameters passed to the script
$projectID = $_GET['project_id'];
$pageNum = $_GET['pageNum'];
$pageSize = 100;
$sort = 'none';

//check if we have abundance for this project
$query = "SELECT has_abundance FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$haveAbundance = false;
$row = mysql_fetch_array($result);
if($row['has_abundance'] == 'Y'){
	$haveAbundance = true;
}

if(isset($_GET['filter'])){
	$filter = $_GET['filter'];	
	$filterValue = $_GET['filterValue'];
}
else{
	$filter = 'none';	
}

if(isset($_GET['sortVal'])){
	$sortVal = $_GET['sortVal'];	
}
else{
	$sortVal = '1';	
}

checkPermission($projectID);

#determine the number of pages

include("../../../includes/slimebase2/browsing/neuroBaseBrowseUtils.php");
displaySeqHeader();
displayAllSeqs();
echo("<div id='bottom-page-nav'>");
displayPageNav("browse");
echo("</div>");

?>