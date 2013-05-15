<?php
include_once('dbConnect.php');
include_once('dbUtils.php');
mysql_select_db('moroz_lab');

//check if the user_id variable is set, if it isn't, assign the visitor to the public user account
if(isset($_SESSION['user_id'])){
	//only check if this is not the public account
	if($_SESSION['user_id'] != 2){
		//make sure the user's session_id matches the id in the sessions table
		$query = "SELECT session_id FROM sessions WHERE user_id ='" . $_SESSION['user_id'] . "'  AND session_id ='" . session_id() . "'";
		$result = mysql_query($query);
		//if it does, update with their latest access
		if(mysql_num_rows($result)){
			$query = "UPDATE sessions SET last_access = CURRENT_TIMESTAMP() WHERE session_id = '" . session_id() . "'";
			mysql_query($query);
		}
		else{
			//using an expired session, regenerate it and insert it into session table
			restartSession();
			$_SESSION['user_id'] = 2;
			$query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '2', '0', CURRENT_TIMESTAMP())";
			$result = mysql_query($query);
		}
	}
	else{
		if(isset($_SESSION['group_user'])){
			//change the user over from public to the user for the group
			$query = "DELETE FROM sessions WHERE session_id='" . session_id() . "'";
			mysql_query($query);
			$_SESSION['user_id'] = $_SESSION['group_user'];
			$query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '" . $_SESSION['group_user'] . "', '0', CURRENT_TIMESTAMP())";
			$result = mysql_query($query);
		}
		else{
		//make sure the session still exists in the db
		$query = "SELECT * FROM sessions WHERE session_id = '" . session_id() . "'";
		$result = mysql_query($query);
		if(mysql_num_rows($result)){
			//update the session's access in the db
			$query = "UPDATE sessions SET last_access = CURRENT_TIMESTAMP() WHERE session_id = '" . session_id() . "'";
			$result = mysql_query($query);		
		}
		else{
			//regenerate session_id and insert the new one in the database
			restartSession();
			$_SESSION['user_id'] = 2;
			$query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '2', '0', CURRENT_TIMESTAMP())";
			$result = mysql_query($query);
		}
		}
		
	}
}
else{
	if(isset($_SESSION['group_user'])){
			$_SESSION['user_id'] = $_SESSION['group_user'];
			$query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '" . $_SESSION['group_user'] . "', '0', CURRENT_TIMESTAMP())";
			$result = mysql_query($query);
	}
	else{
	 //assign this user to the public username and give them a session_id
	 $_SESSION['user_id'] = 2;
	 //add an entry in the database this session
	 $query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '2', '0', CURRENT_TIMESTAMP())";
	 mysql_query($query);
	}
}

if(isset($_SESSION['user_id'])){
	//determine if this user has universal access
	$query = "SELECT user_id FROM permissions WHERE user_id='" . $_SESSION['user_id'] . "' AND project_id =-1";
	$result = mysql_query($query);
	$universalAccess = false;
	while($row = mysql_fetch_array($result)){
		$universalAccess = true;
	}
	$groupMembership = "";
	//check if this user belongs to a group
	$query = "SELECT group_id FROM user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
			$groupMembership .= $row['group_id'];
	}
}
else{
	
	$_SESSION['user_id'] = 2;
	$universalAccess = false;
	$groupMembership = "";
}
//call to gc routine
gc();

function curPageURL() {
 $pageURL = 'http';
 if ($_SERVER["HTTPS"] == "on") {$pageURL .= "s";}
 $pageURL .= "://";
 if ($_SERVER["SERVER_PORT"] != "80") {
  $pageURL .= $_SERVER["SERVER_NAME"].":".$_SERVER["SERVER_PORT"].$_SERVER["REQUEST_URI"];
 } else {
  $pageURL .= $_SERVER["SERVER_NAME"].$_SERVER["REQUEST_URI"];
 }
 return $pageURL;
}

?>