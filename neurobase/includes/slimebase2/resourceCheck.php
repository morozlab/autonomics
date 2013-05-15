<?php

//check if we have homology data for this project
$query = "SELECT sb_id FROM sorted_homology WHERE project_id ='" . $projectID . "' AND sort_id='1' LIMIT 0,1";
$result = mysql_query($query);
$haveAnnotation = false;
if(mysql_num_rows($result)){
	$haveAnnotation = true;
}


//check our annotations
$haveKEGG = false;
$haveGO = false;
$havePfam = false;
$haveAbundance = false;
$query = "SELECT kegg, pfam, go, has_abundance FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$haveKEGG = false;
$row = mysql_fetch_array($result);
if($row['kegg'] == 'Y'){
	$haveKEGG = true;	
}
if($row['pfam'] == 'Y'){
  $havePfam = true;
}

if($row['go'] == 'Y'){
  $haveGO = true;
}
if($row['has_abundance'] == 'Y'){
  $hasAbundance = true;
}

if(isset($sort)){
	if(($haveAbundance == false) && ($sort == 'abundance')){
		$sort = 'evalue';
		$_SESSION['sort'] = 'evalue';
	}
}
?>