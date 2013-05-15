<?php
include("dbConnect.php");
session_start();
mysql_select_db("moroz_lab") or die("Cannot connect to database!");
$query = "SELECT user_id FROM sessions WHERE session_id ='" . session_id() . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$userID = $row['user_id'];
}
if(($userID != $_SESSION['user_id']) || (!(isset($_SESSION['user_id']))) || (session_id() == "")){
	echo("<script language='javascript'>alert('You are currently not logged in, redirecting to log in page!');");
	echo("location.href='http://10.41.128.72/lab_members/member_log_in.html';</script>");	
}
else{
	//update the acess time of the session
	$query = "UPDATE sessions SET last_access = CURRENT_TIMESTAMP() WHERE session_id ='" . session_id() . "'";
	mysql_query($query);
}
?>