<?php

$currentPage = $_SERVER['PHP_SELF'];
$path = explode("/", $currentPage);
$thisPage = $path[sizeof($path) - 1]; 
//check if this was a self-referral
if(!(isset($_SERVER['HTTP_REFERER']))){
	$self = true;
}
else{
	if((preg_match("/\/viewSeqs.php/", $_SERVER['HTTP_REFERER']))){
			$self = true;
	}
	else{
		$self = false;
	}
}

//get the projectID
if(isset($_GET['projectID'])){
	$_SESSION['projectID'] = $_GET['projectID'];
	$projectID = $_SESSION['projectID'];
}
else{
	if(isset($_SESSION['projectID'])){
		$projectID = $_SESSION['projectID'];
	}
	else{
		$projectID = -1;
	}
}

//check if there is a query, i.e. is this a search or not
if(isset($_GET['query'])){
	$search = $_GET['query'];
	if(isset($_SESSION['filter'])){
	$_SESSION['filter'] == 'none';
	}
}

//see if they're trying to search a specific project
if(isset($_GET['searchProj'])){
	$toSearch = $_GET['searchProj'];
	$projectID = $toSearch;
}

if(isset($_SESSION['numAnnots'])){
	$numberAnnotations = $_SESSION['numAnnots'];
}
else{
	$numberAnnotations = 5;
}
//check if the display option for sequence is set
if(isset($_SESSION['display'])){
	$display = $_SESSION['display'];
}
//default display is annotation
else{
	$display = "annotation";
}
//check if filter option is set
if(isset($_GET['filter'])){
	$_SESSION['filter'] = $filter = $_GET['filter'];
	if(isset($_GET['filterValue'])){
		$filterValue = $_GET['filterValue'];
		$_SESSION['filterValue'] = $filterValue;
	}
}
else{
	if(isset($_SESSION['filter'])){
		$filter = $_SESSION['filter'];
		if($filter == "evalue" || $filter == "annotationSource"){
			$filterValue = $_SESSION['filterValue'];
		}
	}
	else{
		$filter = 'none';
	}
	//default visit from other page to not have any filtering
	if($self == false){
		$filter = 'none';
		$_SESSION['filter'] = 'none';
	}
}

$sort = 'evalue';
$sortMod;
if(isset($_GET['sort'])){
	$_SESSION['sort'] = $sort = $_GET['sort'];
}
else{
	if(isset($_SESSION['sort'])){
		$sort = $_SESSION['sort'];
	}
	else{
		$sort = 'evalue';
	}
	if($self == false){
		$sort = 'evalue';
		$_SESSION['sort'] = 'evalue';
	}
}



?>