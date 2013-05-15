<?php
include('../../../includes/dbConnect.php');
include('../../../includes/checkSession.php');
mysql_select_db('moroz_lab');
$projName = $_POST['projName'];
$projParent = $_POST['projParent'];
$query = "SELECT projectID FROM project_directory WHERE project_name = '" . $projName . "' AND child_of ='" . $projParent . "''";
$result = mysql_query($query);
if(mysql_num_rows($result) > 0){
	echo('<script language = "javascript">alert("Project already exists in the selected folder. Please change either the project name or the destination folder.");</script>');
}
else{
	//get the directoy of the parent
	$query = "SELECT path FROM project_directory WHERE projectID = '" . $projParent . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$phpPath = $row['path'] . "\\" . $projName;
		$mysqlPath = addslashes($row['path']) . "\\\\" . $projName;
	}
	mkdir($phpPath) or die("Could not create project folder!");
	//add project to DB
	$query = "INSERT INTO project_directory (project_name, child_of, path, seq_table, num_NT_seqs, num_AA_seqs) VALUES ('" . $projName . "', '" . $projParent . "', '" . $mysqlPath . "', '" . $projName . "_sequences', '0', '0')";
	$result = mysql_query($query);
}

//get path to the newly created project
/*echo('<script language="javascript">alert("' . $path . '");</script>');
//create the folder to hold the project
*/echo('<script language = "javascript">location.href="http://10.41.128.72/admin/seq_view/manage_viewer_db.php";</script>');
?>
