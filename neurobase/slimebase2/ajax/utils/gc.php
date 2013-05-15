<?php
include_once('../../../includes/dbConnect.php');
include('../../../includes/dbUtils.php');
mysql_select_db('moroz_lab');
$selector = "WHERE last_access < date_sub(CURRENT_TIMESTAMP(), interval 1 hour)";
//first, get the list of sessions that are eligible for garbage collection
$query = "SELECT session_id FROM sessions " . $selector;
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	//remove all the important stuff for these sessions
	$session = $row['session_id'];
	$query = "DELETE FROM sessions WHERE session_id ='" . $session . "'";
	mysql_query($query);
	$query = "SELECT blast_id FROM finished_jobs WHERE session_id ='" . $session . "'";
	$result2 = mysql_query($query);
	while($row = mysql_fetch_array($result2)){
		//for each blast job, delete the folder
		removeDir("C:/wamp/www/slimebase2/results/" . $row['blast_id']);
	}	
	$query = "DELETE FROM finished_jobs WHERE session_id = '" . $session . "'";
	mysql_query($query);
	$query = "DELETE FROM active_hits WHERE session_id ='" . $session . "'";
	mysql_query($query);
	$query = "DELETE FROM active_queries WHERE session_id = '" . $session . "'";
	mysql_query($query);
}

?>