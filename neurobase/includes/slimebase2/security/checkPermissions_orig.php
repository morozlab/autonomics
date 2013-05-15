<?php

//check if this use has ability to access the project they're currently viewing
if($universalAccess == false && $projectID != -1){
	$query = "SELECT user_id FROM permissions WHERE project_id ='" . $projectID . "' and user_id = '" . $_SESSION['user_id'] . "'";
	$result = mysql_query($query);
	if(!(mysql_num_rows($result))){
		//check if there is a group associated to this user ID that has access to this project
		$query = "SELECT user_id FROM user_group_membership, project_group_membership WHERE user_group_membership.group_id = project_group_membership.group_id AND user_group_membership.user_id='" . $_SESSION['user_id'] . "' AND project_group_membership.project_id='" . $projectID . "'";
		$result2 = mysql_query($query);
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
				$rest = "slimebase2/login.php?msg=permfail&projectID=" . $projectID;
			}
			header("Location: " . $base . $rest);
	/*echo('<script language="javascript">location.href="http://74.252.103.104:8888/slimebase2/login.php?msg=permfail&projectID=' . $projectID . '";</script>');
*/		}
	}
}

function checkPermission(){
	global $unviersalAccess;
	global $projectID; 
	if($universalAccess == false && $projectID != -1){
		$query = "SELECT user_id FROM permissions WHERE project_id ='" . $projectID . "' and user_id = '" . $_SESSION['user_id'] . "'";
		$result = mysql_query($query);
		if(!(mysql_num_rows($result))){
			//check if there is a group associated to this user ID that has access to this project
			$query = "SELECT user_id FROM user_group_membership, project_group_membership WHERE user_group_membership.group_id = project_group_membership.group_id AND user_group_membership.user_id='" . $_SESSION['user_id'] . "' AND project_group_membership.project_id='" . $projectID . "'";
			$result2 = mysql_query($query);
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
					$rest = "slimebase2/login.php?msg=permfail&projectID=" . $projectID;
				}
				header("Location: " . $base . $rest);
	/*echo('<script language="javascript">location.href="http://74.252.103.104:8888/slimebase2/login.php?msg=permfail&projectID=' . $projectID . '";</script>');
*/			}
		}
	}
}


?>