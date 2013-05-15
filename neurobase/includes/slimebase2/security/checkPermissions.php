<?php

$page = $_SERVER['PHP_SELF'];
$array = explode("/", $page);
$pageName = $array[count($array) - 1];

$projectID = $_GET['projectID'];
checkPermission();


function checkPermission(){
	global $universalAccess;
	global $projectID; 
	if($universalAccess == false){
		$query = "SELECT permissions.user_id FROM permissions JOIN sessions ON permissions.user_id=sessions.user_id WHERE permissions.project_id ='" . $projectID . "' and session_id = '" . session_id() . "'";
		$result = mysql_query($query);
		if(!(mysql_num_rows($result))){
			//check if there is a group associated to this user ID that has access to this project
			$query2 = "SELECT user_id FROM user_group_membership, project_group_membership WHERE user_group_membership.group_id = project_group_membership.group_id AND user_group_membership.user_id='" . $_SESSION['user_id'] . "' AND project_group_membership.project_id='" . $projectID . "'";
			$result2 = mysql_query($query2);
			if(!(mysql_num_rows($result2))){
			
				if(isset($_SERVER['HTTP_HOST'])){
					$base = "http://" . $_SERVER['HTTP_HOST'] . "/";
				}
				else{
					$base = "http://74.252.103.104:8888/";
				}
				if($_SESSION['user_id'] == 4){
					//special case for reviewers
					$rest = "slimebase2/review_login.php?msg=permfail&projectID=" . $projectID;
				}
				else{
					//everyone else
					$rest = "slimebase2/login.php?msg=permfail&projectID=" . $projectID  . "&check=" . $query;
				}
				header("Location: " . $base . $rest);
	/*echo('<script language="javascript">location.href="http://74.252.103.104:8888/slimebase2/login.php?msg=permfail&projectID=' . $projectID . '";</script>');
*/			}
		}
	}
}


?>