<?php

session_start();

include('../../includes/dbConnect.php');
mysql_select_db('moroz_lab');

$program = $_GET['program'];

//check for universal access
$query = "SELECT user_id FROM permissions WHERE user_id='" . $_SESSION['user_id'] . "' AND project_id =-1";
$result = mysql_query($query);
$universalAccess = false;
while($row = mysql_fetch_array($result)){
	$universalAccess = true;
}

$groupMemberhsip = "";
//check if this user belongs to a group
$query = "SELECT group_id FROM user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$groupMembership .= $row['group_id'];
}

$toEcho = "";
if($program == "blastn" || $program == "tblastn" || $program == "tblastx"){
	if($universalAccess == true){
		$query = "SELECT projectID, project_name FROM project_directory WHERE num_NT_seqs != '0' ORDER BY project_name";
	}
	else{
		if($groupMembership != ""){
			$query = "SELECT projectID, project_name FROM project_directory, project_group_membership, user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "' AND user_group_membership.group_id = project_group_membership.group_id AND project_directory.projectID = project_group_membership.project_id";
		}
		else{
			$query = "SELECT projectID, project_name FROM project_directory, users, permissions WHERE num_NT_seqs != '0' AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
		}
	}
}
else{
	if($universalAccess == true){
		$query = "SELECT projectID, project_name FROM project_directory WHERE num_AA_seqs != '0' ORDER BY project_name";
	}
	else{
		$query = "SELECT projectID, project_name FROM project_directory, users, permissions WHERE num_AA_seqs != '0' AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
	}
}
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$toEcho = $toEcho . "<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>"; 
}
echo($toEcho);

?>