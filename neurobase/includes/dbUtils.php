<?php


function restartSession(){
	$_SESSION = array();
	if (isset($_COOKIE[session_name()])) {
    	setcookie(session_name(), '', time()-42000, '/');
	}
	session_regenerate_id();
}

//checks the sessions table and removes sessions that are eligible for collection (currently three days old
function gc(){
	$selector = "WHERE last_access < date_sub(CURRENT_TIMESTAMP(), interval 1 hour)";
	//first, get the list of sessions that are eligible for garbage collection
	$query = "SELECT session_id FROM sessions " . $selector;
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		//remove all the important stuff for these sessions
		$session = $row['session_id'];
		$query = "DELETE FROM sessions WHERE session_id ='" . $session . "'";
		mysql_query($query);
		$query = "DELETE FROM active_hits WHERE session_id ='" . $session . "'";
		mysql_query($query);
		$query = "DELETE FROM active_queries WHERE session_id = '" . $session . "'";
		mysql_query($query);
		$query = "DELETE FROM sequence_cart WHERE session_id ='" . $session . "'";
		mysql_query($query);
		$query = "SELECT blast_id FROM finished_jobs WHERE session_id ='" . $session . "'";
		$result2 = mysql_query($query);
		while($row = mysql_fetch_array($result2)){
			//for each blast job, delete the folder
			removeDir("/var/www/html/neurobase/slimebase2/results/" . $row['blast_id']);
		
		}
		$query = "DELETE FROM finished_jobs WHERE session_id = '" . $session . "'";
		mysql_query($query);	
	}
}

function removeDir($dir)
{
	if (!file_exists($dir)) return true;
	if (!is_dir($dir) || is_link($dir)) return unlink($dir);

	foreach (scandir($dir) as $item)
	{
		if ($item == "." || $item == "..") continue;
		if (!removeDir($dir . "/" . $item))
		{

			if (!removeDir($dir . "/" . $item)) return false;
		};
	}
	return rmdir($dir);
}


?>