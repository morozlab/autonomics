<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<?php 
include('../includes/dbConnect.php');
mysql_select_db('moroz_lab');
if(isset($_GET['projectID'])){
	$projectID = $_GET['projectID'];
}
else{
	$projectID = 57;
}
if(isset($_GET['query'])){
	$search = $_GET['query'];
}
else{
	//get the project name
	$query = "SELECT project_name FROM project_directory WHERE projectID='" . $projectID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$projectName = $row['project_name'];
}
if(isset($_GET['pageSize'])){
	$pageSize = $_GET['pageSize'];
}
else{
	$pageSize = 20;
}
if(isset($_GET['curPage'])){
	$curPage = $_GET['curPage'];
}
else{
	$curPage = 1;
}
 

?>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php 
if(isset($projectName)){
	echo("Sequences - " . $projectName);
}
else{
	echo("Search results - " . $search);
}
?> </title>
<link href="../css/slimebase.css" rel="stylesheet" type="text/css" />
</head>

<body>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<?php include('../includes/slimebaseNav.php'); ?>
        
    <div id="search">
        	<form name='search' id='searchForm' method='post' action='searchResults.php'>
            	<input type='text' size='20' name='query' value='search...' onFocus='this.select();' />
                <input type='hidden' name='ref' value='browse.php' />
            	<select name='database' style='width: 150px;'>
                <option value='-1' selected="selected">Select Database</option>
				<?php
					$query = "SELECT project_name, projectID FROM project_directory WHERE (num_AA_seqs != 0) OR (num_NT_seqs != 0)";
					$result = mysql_query($query);
					while($row = mysql_fetch_array($result)){
						echo("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
					}
				?>
                </select>
            </form>
        </div>
    </div>
    <div id='mainContent'>
    	<?php
			//first case, there was no search, so just display all of the sequences for the project
			if(!(isset($search))){
				displayAllSeqs($projectID, $pageSize, $curPage);				
			}	
		?>
    </div>   
    </div>
</div>
</body>
</html>

<?php

function displayAllSeqs($projectID, $pageSize, $curPage){
	$query = "SELECT num_AA_seqs, num_NT_seqs, project_name FROM project_directory WHERE projectID ='" . $projectID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	echo("<div id='resultHeader'>" . $row['project_name'] . "</div>");
	if($row['num_NT_seqs'] > $row['num_AA_seqs']){
		$numSeqs = $row['num_NT_seqs'];
	}
	else{
		$numSeqs = $row['num_AA_seqs'];
	}
	if($numSeqs == 0){
		
	}
	else{
		//determine the pagination
		$last = ceil($numSeqs/$pageSize);
		if($curPage < 1){
			$curPage = 1;
		}
		else if($curPage > $last){
			$curPage = $last;
		}
		$pageOffset = ($curPage - 1) * $pageSize;
		$limiter = "LIMIT " . $pageOffset . ", " . $pageSize;
		$seqCounter = $pageOffset + 1;
		$query = "SELECT sb_id, description, length, date FROM " . $projectID . "_sequences " . $limiter;
		$result = mysql_query($query);
		//display page statistics and navigations
		echo("<div id='pageInfoContainer'>");
		echo("<div id='resultDescription'>");
		echo("Sequences <b>" . ($pageOffset  + 1) . "</b> - <b>" . ($pageOffset + $pageSize) . "</b> of " . $numSeqs);
		echo("</div>");
		echo("<div id='resultNav'>");
		if($curPage != $last){
			echo("<div id='resultLink'><a href='viewSeqs.php?projectID=" . $projectID . "&pageSize=" . $pageSize . "&curPage=" . ($curPage + 1) . "'>next</a></div>");
		}
		//show the page input box
		echo("<div id='resultForm'><form name='pageSelect' action='viewSeqs.php' method='get'>");
		echo("page <input type='text' size='" . strlen($curPage) . "' name='curPage' value='" . $curPage . "' onFocus='this.select();'/> of " . $last);
		echo("<input type='hidden' name='pageSize' value='" . $pageSize . "'/>");
		echo("<input type='hidden' name='projectID' value='" . $projectID . "' />");
		echo("</form></div>"); 
		if($curPage != 1){
			echo("<div id='resultLink'><a href='viewSeqs.php?projectID=" . $projectID . "&pageSize=" . $pageSize . "&curPage=" . ($curPage - 1) . "'>previous</a></div>");	
		}
		echo("</div>");
		echo("</div>");
		//iterate over the sequences and display the information for each of them
		while($row = mysql_fetch_array($result)){
			echo("<div id='itemEntry'>");
			echo("<div id='seqNumber'><b>" . $seqCounter . ".</b></div>");
			echo("<div id='seqDescription'>");
			//split the sb out from the 
			$description = $row['description'];
			$splitDesc = split("\|", $description);
			$sb = $splitDesc[1];
			$description = $splitDesc[2];
			echo("sb|<a href='seqDetail.php?sbid=" . $sb . "'>" . $sb . "</a>|" . $description);
			echo("</div>");
			//reformat the date from the db
			$splitDate = split("\-", $row['date']);
			$date = $splitDate[1] . "/" . $splitDate[2] . "/" . $splitDate[0];
			echo("<div id='seqDate'>");
			echo("modified: " . $date);
			echo("</div>");
			echo("</div>");
			$seqCounter++;
		}
	}
}
