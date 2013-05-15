<?php 
session_start();
include('../../includes/dbConnect.php');
mysql_select_db('moroz_lab');
if(isset($_GET['projectID'])){
	$projectID = $_GET['projectID'];
}
else{
	$projectID = -1;
}

include('../../includes/restrictAccess.php');
include('../../includes/slimebase2/security/checkPermissions.php');

//get the project name
$query = "SELECT project_name FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$projectName = $row['project_name'];


?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php echo("GOBrowser - " . $projectName); ?></title>
<link href="../../css/slimebase_basic.css" rel="stylesheet" type="text/css" />
<link href="../../css/slimebaseGO.css" rel="stylesheet" type="text/css" />
</head>

<body>
<script type="text/javascript" src="../javascript/mochikit/lib/MochiKit/MochiKit.js"></script>
<script type="text/javascript" src="../javascript/plotkit/PlotKit/Base.js"></script>
<script type="text/javascript" src="../javascript/plotkit/PlotKit/Layout.js"></script>
<script type="text/javascript" src="../javascript/plotkit/PlotKit/Canvas.js"></script>
<script type="text/javascript" src="../javascript/plotkit/PlotKit/SweetCanvas.js"></script>
<script type="text/javascript" src="../javascript/plotkit/PlotKit/excanvas.js"></script> 
<script type='text/javascript' src="../javascript/goGraphs.js"></script>

<script language='javascript' type='text/javascript'>
	var projectID = <?php echo($projectID); ?>;


</script>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<img style='float:left; margin-right: 20px;' src="../images/slimeicon.jpg" alt="SlimeBase" />
        <div class='navElement'><a href='../main.php'>Home</a></div>
        <div class='navElement'><?php
			if(isset($groupMembership) && ($groupMembership != "")){
				if($groupMembership == 2){
					?> <a href='../nar.php'>Projects</a> <?php	
				}
				else if($groupMembership == 1){
					?> <a href='../chpnas.php'>Projects</a><?php	
				}
			}
			else{ ?>
        		<a href='../browse.php'>Projects</a>
        	<?php } ?></div>
        <div class='navElement'><a href='../blast.php'>BLAST</a></div>
        
    <div id="search">
        	<form name='search' id='searchForm' method='get' action='../viewSeqs.php'>
            	search <input type='hidden' name='ref' value='browse' /> 
                <select name='projectID' style='width: 150px;'>
                <option value='-1' selected="selected">Select Database</option>
				<?php
				if($universalAccess == false){
					if($groupMembership != ""){
						$query = "SELECT projectID, project_name FROM project_directory, project_group_membership, user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "' AND user_group_membership.group_id = project_group_membership.group_id AND project_directory.projectID = project_group_membership.project_id";
					}
					else{
						$query = "SELECT projectID, project_name, user_name FROM project_directory, permissions, users WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
					}
				}
				else{
					$query = "SELECT projectID, project_name FROM project_directory WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) ORDER BY project_name";
				}
					$result = mysql_query($query);
					while($row = mysql_fetch_array($result)){
						echo("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
					}
				?>
                </select>
                for <input type='text' size='20' name='query' /> <input type='submit' name='submit' value='Go' />
</form>
        </div>
    </div>
    <?php include('../../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='goNav'><h1>NeuroBase GO Browser</h1></div>
        <div id='goHeader'></div>
        <div id='goContent'>
        	<div id='goDescription'>GO Categories - <?php echo($projectName); ?></div>
            <div id='top_level'>
        	</div>
            <div id='summaryDivider'><h1>GO Term Distribution</h1></div>
			<!--
            <div id='goDisplayOptions'>
            	<table width="898">
                	<tr><td>display:</td><td><input type='checkbox' checked="checked" name='top_level' onclick="toggleDisplay(this.name);" />top-level GO categories</td><td><input type='checkbox' name='bio_p' onclick='toggleDisplay(this.name);' />biological process</td><td><input type='checkbox' name='mollec' onclick="toggleDisplay(this.name);" />mollecular function</td><td><input type='checkbox' name='cell' onclick='toggleDisplay(this.name);' />cellular component</td></tr>
                </table>
            </div>
            -->
            <div id='goDisplayOptions'><div class='option' style='padding-top: 3px;'>display:</div><div class='option'><input type='checkbox' checked="checked" id='top_level_check' name='top_level' onclick="toggleDisplay(this.name);" />top-level GO categories</div>
        
            <div class='option'><input id='bio_p_check' type='checkbox' name='bio_p' value='checked' onclick='toggleDisplay(this.name);' />biological process</div>
           	 <div class='option'><input id='mollec_check' type='checkbox' name='mollec' onclick="toggleDisplay(this.name);" />molecular function</div>
           	 <div class='option'><input id='cell_check' type='checkbox' name='cell' onclick='toggleDisplay(this.name);' />cellular component</div>
            </div>
            
            <div id='goGraphs'>
            	<h1 id='top_level_h1' style="float:left; width:700px; margin-top:40px; margin-bottom: 30px; padding-left: 20px; font-family: Geneva, Arial, Helvetica, sans-serif;	font-size: 1.2em;color:#2F738A;">Top-Level GO Categories</h1>
            <div id='top_level_div' style="float:left; width: 300px; margin-left: 30px;"><canvas id="top_categories_canvas" height="300" width="300"><div style="width: 500px;"></div></canvas></div>
        	</div>
            <h1 id='bio_p_h1' class='goTitle'>Biological Process</h1>
            <div id='bio_p_div' class="graphContainer"><canvas id="bio_p_canvas" height="250" width="800"></canvas></div>
				<h1 id='mollec_h1' class='goTitle'>Molecular Function</h1>
            <div id='mollec_div' class="graphContainer"><canvas id="mollec_canvas" height="250" width="800"></canvas></div>
				<h1 id='cell_h1' class='goTitle'>Cellular Component</h1>
            <div id='cell_div' class="graphContainer"><canvas id="cell_canvas" height="250" width="800"></canvas></div>
        </div>
        <div id='goFooter'></div>
        <script language='javascript' type='text/javascript'>
			//initialize the top level GO categories, if we're viewing a project
			if(projectID != -1){
				initTopLevel();
				MochiKit.DOM.addLoadEvent(gatherGraphData);
			}
		</script>
    </div>   
    </div>
</div>
</body>
</html>
