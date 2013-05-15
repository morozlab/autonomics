<?php
include('../../../includes/dbConnect.php');
include('../../../includes/checkSession.php');
mysql_select_db('moroz_lab');
$projID = $_POST['projID'];
//remove remote folder for project
$query = "SELECT path FROM project_directory WHERE projectID = '" . $projID . "'";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$path = $row['path'];
}
if(!(preg_match('!C:\\\\wamp!', $path))){
	die("SYSTEM ERROR: Stopped potential fatal delete outside of WAMP directory!");
}
$command = "RD /S /Q " . $path;
system($command);
//remove the entry in the project directory
$query = "DELETE FROM project_directory WHERE projectID = '" . $projID . "'";
$result = mysql_query($query);
//remove entries from the project file table
$query = "DELETE FROM project_files WHERE projectID = '" . $projID . "'";
$result = mysql_query($query);
//remove sequence tables for the project
$query = "DROP TABLE " . $projID . "_sequences";
$result = mysql_query($query);
echo('<script language = "javascript">location.href="http://150.176.130.196:8888/admin/seq_view/manage_viewer_db.php";</script>');

?>