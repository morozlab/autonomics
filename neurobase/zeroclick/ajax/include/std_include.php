<?php

session_start();
$_GET['projectID'] = -123109123;
$path = $_SERVER['DOCUMENT_ROOT'];
include_once($path . '/includes/restrictAccess.php');
include_once($path . '/includes/slimebase2/security/checkPermissions.php');
include_once($path . '/zeroclick/ajax/db.php');
include_once($path . '/zeroclick/ajax/iondb.php');


function sanitize_checkbox($checked){
	if($checked == "true"){
		return true;
	}
	else{	
		return false;
	}
		
}

function sanitize_number($num){
	
	if(is_numeric($num)){
		return $num;	
	}
	else{
		return null;	
	}
	
}

?>