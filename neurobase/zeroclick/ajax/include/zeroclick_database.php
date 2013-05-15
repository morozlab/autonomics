<?php

include_once('std_include.php');

function job_exists($pid, $job_type){
	$dbh = make_db_connection();
	$sth = $dbh->prepare("SELECT * FROM jn_mapping WHERE project_id= ? AND job_type = ?");
	$sth->execute(array($pid, $job_type));
	$jid = '';
	while($row = $sth->fetch()){
		$jid = $row['job_id'];	
	}
	if($jid == ''){
		return false;		
	}
	return true;	
}

function job_configured($pid, $job_type){
	$dbh = make_db_connection();
	$sth = $dbh->prepare("SELECT * FROM configuration WHERE project_id=? AND job_type=?");
	$sth->execute(array($pid, $job_type));
	while($row = $sth->fetch()){
		return true;	
	}
	return false;
}

?>