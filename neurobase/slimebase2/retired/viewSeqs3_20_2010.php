<?php 
session_start();
include('../includes/dbConnect.php');
include('../includes/slimebase2/seqViewOptionCheck.php');
include('../includes/slimebase2/seqViewResourceCheck.php');
include('../includes/restrictAccess.php');
include('../includes/slimebase2/security/checkPermissions.php');
mysql_select_db('moroz_lab');

//get the project name
if($projectID != -1){
	$query = "SELECT project_name FROM project_directory WHERE projectID='" . $projectID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$projectName = $row['project_name'];
}
else{
	$projectName = "SlimeBase Search";
}
if(isset($_GET['pageSize'])){
	$pageSize = $_GET['pageSize'];
}
else{
	$pageSize = 100;
}

include('../includes/slimebase2/browsing/browseUtilsv2.php');

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php if(isset($projectName)){	echo("Sequences - " . $projectName);}else if(isset($search)){	echo("Search results - " . $search);}else{ echo("SlimeBase Search");}?></title>

<link href="../css/seq_view.css" rel="stylesheet" type="text/css" />
<link href="../css/slimebase_white.css" rel="stylesheet" type="text/css" />
</head>

<body>
<script language='javascript' type='text/javascript'>
sfHover = function() {
	var sfEls = document.getElementById("projectOptions").getElementsByTagName("LI");
	for (var i=0; i<sfEls.length; i++) {
		sfEls[i].onmouseover=function() {
			this.className+=" sfhover";
		}
		sfEls[i].onmouseout=function() {
			this.className=this.className.replace(new RegExp(" sfhover\\b"), "");
		}
	}
}
if (window.attachEvent) window.attachEvent("onload", sfHover);
</script>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<?php include('../includes/slimebaseNav.php'); ?>
        
    <div id="search">
        	<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        </div>
    </div>
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<?php displayProjectOptions(); ?>
    	<?php
			//first case, the user is wanting to search a specific project
			if(isset($toSearch)){
				displaySearchOptions($toSearch);
			}
			//second case, we're viewing sequences from a search result
			else if(isset($search)){
				if($projectID == -1){
					displayError('noDB');
				}
				else{
					displaySearchResults();
				}
			}
			//third case, user is browsing entire project
			else{
				displayAllSeqs();
				displayPagination();
			}
		?>
    </div>   
    </div>
</div>
</body>
</html>


