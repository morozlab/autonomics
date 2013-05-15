<?php 
session_start();
$_SESSION['group_user'] = 13;
include('../includes/dbConnect.php');
include('../includes/restrictAccess.php');
include('../includes/slimebase2/projects/checkProjectResources.php');
mysql_select_db('moroz_lab');
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>NeuroBase: A Comparative Neurogenomics Database</title>
<link href="../css/neurobase_basic.css" rel="stylesheet" type="text/css" />
<link href="../css/browse_new.css" rel="stylesheet" type="text/css" />
<script language='javascript' type='text/javascript' src='javascript/mootools/mootools.js'></script>
<script language='javascript' type='text/javascript' src='javascript/projects/projectUtils.js'></script>
<script language='javascript' type='text/javascript' src="javascript/seqUtils.js"></script>
<script language='javascript' type='text/javascript'>
	var detailExpanded = false;
	var currentSelected = "";
</script>
</head>

<body onload='positionHeader();'>
<div id='neurobase-logo'><img src="images/neuro_base_logo.png" width="624" height="69" /></div>
        <div id='neurobase-nav-left'><img src="images/neuro_base_nav_left.png" /></div>
        <div id='neurobase-nav'><div class='first-nav'><a href='../index.php'>home</a></div><div class='nav-element'><?php
			if(isset($groupMembership) && ($groupMembership != "")){
				if($groupMembership == 2){
					?> <a href='nar.php'>projects</a> <?php	
				}
				else if($groupMembership == 1){
					?> <a href='chpnas.php'>projects</a><?php	
				}
			}
			else{ ?>
        		<a href='browse.php'>projects</a>
        	<?php } ?></div><div class='nav-element'><a href='blast.php?view=blast'>BLAST</a></div><div class='nav-element'><a href='viewSeqs.php?searchProj=-1'>search</a></div></div>
<div id='container'>
	<div id="content">
    <?php //include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='selectContainer'>
            <div id='selectBoxMiddle'>
            	<div id='selectBox'>
                	<form id='projectSelect' name="projectSelect">
                    	<select id='projectSelectElement' name='projectSelectElement' onchange="loadProjectDetails();">
                        	<option value='-1'>Select species & tissue here</option>
                            <?php
								$query = "SELECT projectID, project_name FROM project_directory, project_group_membership WHERE project_directory.projectID = project_group_membership.project_id AND project_group_membership.group_id='2' ORDER BY project_name";
								$result = mysql_query($query);
								while($row = mysql_fetch_array($result)){
									print("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
								} 
							?>
                        </select>
                    </form> 
                </div>
            </div>
        </div>
        <div id='projDetailContainer'>
        	<div id='projDetailContent'>
            	<div id='organism-description'></div>
            	<div id='projectImage'></div>
                <div id='projectStats'>
                    <div id='pSContent'><h1>Assembly Statistics</h1>
                    	<p id='pSContentParagraph'>
                        sequenced using 454<br />
                        last updated on 6/20/2010<br />
                        7658 contigs<br />
                        765 singlets<br />
                        303bp average contig length<br />
                        6000 annotated sequences<br />
                        </p>
                    </div>
                </div>
            </div>
        </div>
   </div>
   </div> 
</div>
<div id='extraContent'>
	<div id='extraContentNav'>
    </div>
    <div id='extraContentMain'>
    	
    </div>
</div>
<?php
if(isset($_SESSION['project_loaded'])){	
?>
<script type='text/javascript' language="javascript">
	var projSelect = document.getElementById('projectSelectElement');
	projSelect.selectedIndex = 0;
</script>
<?php
}
?>
</body>
</html>
