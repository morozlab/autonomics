<?php

$univeralAccess = false;
$query = "SELECT user_id FROM permissions WHERE project_id='-1' AND user_id='" . $_SESSION['user_id'] . "'";
$result = mysql_query($query);
if(mysql_num_rows($result)){
	$universalAccess = true;	
}

function checkPermission($projectID){
	global $universalAccess;
	$userID = $_SESSION['user_id'];
	if($universalAccess == false && $projectID != -1){
		$query = "SELECT user_id FROM permissions WHERE project_id ='" . $projectID . "' and user_id = '" . $_SESSION['user_id'] . "'";
		$result = mysql_query($query);
		if(!(mysql_num_rows($result))){
			//check if there is a group associated to this user ID that has access to this project
			$query = "SELECT user_id FROM user_group_membership, project_group_membership WHERE user_group_membership.group_id = project_group_membership.group_id AND user_group_membership.user_id='" . $_SESSION['user_id'] . "' AND project_group_membership.project_id='" . $projectID . "'";
			$result2 = mysql_query($query);
			if(!(mysql_num_rows($result2))){
				die;
		/*echo('<script language="javascript">location.href="http://74.252.103.104:8888/slimebase2/login.php?msg=permfail&projectID=' . $projectID . '";</script>');
*/			}
		}		
	}
	//permissions are ok, check if the session is valid
	$query = "SELECT user_id FROM sessions WHERE session_id='" . session_id() . "'";
	$result = mysql_query($query);
	if(!(mysql_num_rows($result))){
		die;
	}
}

?>