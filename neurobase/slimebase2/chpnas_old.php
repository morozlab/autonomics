<?php 
session_start();
$_SESSION['group_user'] = 14;
include('../includes/dbConnect.php');
include('../includes/restrictAccess.php');
include('../includes/slimebase2/projects/checkProjectResources.php');
mysql_select_db('moroz_lab');
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Cancer and Homarus PNAS Comparative Online Database</title>
<link href="../css/slimebase_basic_new.css" rel="stylesheet" type="text/css" />
<link href="../css/browse.css" rel="stylesheet" type="text/css" />
<script language='javascript' type='text/javascript' src='javascript/mootools/mootools.js'></script>
<script language='javascript' type='text/javascript' src='javascript/projects/projectUtils.js'></script>
<script language='javascript' type='text/javascript' src="javascript/seqUtils.js"></script>
<script language='javascript' type='text/javascript'>
	var detailExpanded = false;
	var currentSelected = "";
</script>
</head>

<body>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<?php include('../includes/slimebaseNav.php'); ?>
        
    <div id="search">
        	<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        </div>
    </div>
    <?php //include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='selectContainer'>
        	<div id='selectBoxTextTop'>Getting started is easy.</div>
            <div id='selectBoxMiddle'>
            	<div id='selectBox'>
                	<form id='projectSelect' name="projectSelect">
                    	<select id='projectSelectElement' name='projectSelectElement' onchange="loadProjectDetails();">
                        	<option value='-1'>Select species & tissue here</option>
                            <?php
								$query = "SELECT projectID, project_name FROM project_directory, project_group_membership WHERE project_directory.projectID = project_group_membership.project_id AND project_group_membership.group_id='1' ORDER BY project_name";
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
            	<div id='projectImage'></div>
                <div id='projectStats'>
                	<div class='pSCorner'><img src='images/browse/detail_box_top_left.png' /></div>
                    <div class='pSMiddle'></div>
                    <div class='pSCorner'><img src='images/browse/detail_box_top_right.png' /></div>
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
<?php
}
?>
</body>
</html>
