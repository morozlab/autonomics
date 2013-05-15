<?php

include('../../includes/dbConnect.php');

session_start();

mysql_select_db('moroz_lab');

$session = session_id();

$query = "SELECT session_id FROM sessions WHERE session_id ='" . $session . "'";
$result = mysql_query($query);
if(mysql_num_rows($result)){
	//update the entry in the sessions table
	$query = "UPDATE sessions SET last_access = CURRENT_TIMESTAMP() WHERE session_id = '" . session_id() . "'";
	mysql_query($query);
	echo('yes');
}
else{
	echo('no');
}

?>