<?php
include_once("../includes/dbConnect.php");
include_once("../includes/dbUtils.php");
$username = $_POST['userName'];
$password = $_POST['passWord'];
if(isset($_POST['ref'])){
	$ref = $_POST['ref'];
}
else{
	$ref = 'none';
}
$encryptedPW = md5($password);
$fetchedPW;
$databaseConnection = mysql_select_db("moroz_lab");
if($databaseConnection == false){
	die("Could not connect to database!<br>");
}
$query = "SELECT password, user_id FROM users WHERE user_name='" . $username . "'";
$result = mysql_query($query);
$fetchedPW = 'klojawklijasijeaselkjaseas';
while($row = mysql_fetch_array($result)){
	$fetchedPW = $row['password'];
	$userID = $row['user_id'];
}
#password matches, start session with userID 
if($encryptedPW == $fetchedPW){
	session_start();
	//if this condition is true, the user must have already been logged in...destroy the session and start a new one
	if(isset($_SESSION['user_id'])){
		restartSession();
	}
	else{
	}
	//log the user in
	$query = "INSERT INTO sessions (session_id, user_id, active_job, last_access) VALUES ('" . session_id() . "', '" . $userID . "', '0', CURRENT_TIMESTAMP())";
	$result = mysql_query($query);
	//do some garbage collection
	gc();
	//set session variables for the userID and username
	$_SESSION['user_id'] = $userID;
	$_SESSION['username'] = $username;
	$redirect;
	if($ref == 'none'){
		if($username == "admin"){
			$redirect = '../admin/admin_home.php';
		}
		else{
			$redirect = 'member_homepage.php';
		}
	}
	else{
		$redirect = $ref;
	}
	echo("<script language='javascript'>location.href='" . $redirect . "';</script>");
}
else{
	if($ref == 'none'){
		header("Location: member_log_in.html");
	}
	else{
		if(isset($_SERVER['HTTP_REFERER'])){
			if(preg_match("/\.php?/", $_SERVER['HTTP_REFERER'])){
				header("Location: " . $_SERVER['HTTP_REFERER'] . "&origref=" . $ref . "&msg=credentials");
			}
			else{
				header("Location: " . $_SERVER['HTTP_REFERER'] . "?origref=" . $ref . "&msg=credentials");
			}
		}
		else{
			header("Location: ../slimebase2/login.php?origref=" . $ref . "&msg=credentials");
		}
	}
}
?>


