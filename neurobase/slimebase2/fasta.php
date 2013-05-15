<?php 
session_start();
include('C:/wamp/www/includes/dbConnect.php');
include('C:/wamp/www/includes/restrictAccess.php');
mysql_select_db('moroz_lab');

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php if(isset($projectName)){	echo("Sequences - " . $projectName);}else if(isset($search)){	echo("Search results - " . $search);}else{ echo("SlimeBase Search");}?></title>
<link href="../css/fasta.css" rel='stylesheet' type='text/css' />
<link href="../css/slimebase_white.css" rel="stylesheet" type="text/css" />
<script language='javascript' type='text/javascript' src='javascript/project_view/projViewUtils.js'></script>
<script language='javascript' type='text/javascript' src='javascript/seqUtils.js'></script>
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
    	<?php include('../includes/slimebaseNav_new.php'); ?>
        
    <div id="search">
        	<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        </div>
    </div>
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='fastaHeader'>Custom FASTA</div>
        <div class='fastaOptions'></div>
        <div id='fastaContents'>
        	<table>
            	<?php
					$query = "SELECT * FROM sequence_cart WHERE session_id='" . session_id(). "'";
					$result = mysql_query($query);
					while($row = mysql_fetch_array($result)){
						
					}
				?>
            </table>
        </div>
        <div class='fastaOptions'></div>
        
    </div>   
    </div>
</div>
</body>
</html>


