<?php
session_start();
include('../includes/dbConnect.php');
include('../amfphp/services/Util.php');
include('../includes/restrictAccess.php');
mysql_select_db('moroz_lab');
if(isset($_GET['projectID'])){
	$projectID = $_GET['projectID'];
}
else{
	$projectID = -1;
}
if(isset($_GET['pathwayID'])){
	$pathwayID = $_GET['pathwayID'];	
}
else{
	$pathwayID = -1;
}
include('../includes/slimebase2/security/checkPermissions.php');
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<link rel="icon" type="image/png" href="../favicon.png" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>KEGG Pathway Viewer</title>
<link href="../css/kegg.css" rel="stylesheet" type="text/css" />
<link href='../css/slimebase_basic.css' rel='stylesheet' type='text/css' />
</head>

<body>
<div id='container'>
	<div id="content">
		<div id='header'>
    		<?php include('../includes/slimeBaseNav.php'); ?>        
    		<div id="search">
        		<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        	</div>
    	</div>
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
<!-- load neurobase utility scripts -->
<script type='text/javascript' src="javascript/seqUtils.js"></script>
<!-- load the kegg drawing utilities -->
<script language='javascript' type='text/javascript' src='javascript/ocanvas/ocanvas-2.1.1.js'></script>
<script type='text/javascript' src='javascript/kegg/keggUtils.js'></script>
		<div id='kegg-header'>KEGG Browser</div>
		<div id='kegg-graphics-container'>
			<canvas id='kegg-graphics' class='kegg-graphics'></canvas>
		</div>
        <div id='kegg-pathway-select'><span>Components of</span> <form id='keggSelect' name="keggSelect" action='kegg.php' method='get'>
        				<input type='hidden' name='projectID' value='<?php echo($projectID); ?>'  />
                    	<select id='keggSelectElement' name='pathwayID' onchange="document.keggSelect.submit()">
                            <?php
								$query = "SELECT DISTINCT path_id, kegg_pathway_description FROM kegg_annotations WHERE project_id='" . $projectID . "' ORDER BY kegg_pathway_description";
								$result = mysql_query($query);
								while($row = mysql_fetch_array($result)){
									if($row['path_id'] == $pathwayID){
										print("<option value='" . $row['path_id'] . "' selected='selected'>" . $row['kegg_pathway_description'] . "</option>");
									}
									else{
										print("<option value='" . $row['path_id'] . "'>" . $row['kegg_pathway_description'] . "</option>");
									}

								} 
        //get the sequence type for this project
        $query = "SELECT default_type FROM project_directory WHERE projectID=" . $projectID;
        $result = mysql_query($query);
        while($row = mysql_fetch_array($result)){
          $seq_type = $row['default_type'];
        }

	       			?>
                        </select>
                    </form> 
                    <?php echo("<div id='kegg-detail-options'><div class='kegg-option'><a href='download/doDownload.php?method=exportKEGG&args=" . $projectID . "," . $pathwayID . "," . $seq_type . "'>export sequences</a></div></div>"); ?></div>
        </div>
        <div id='kegg-pathway-components'></div>
	</div>
</div>
<script type='text/javascript'>
	drawKeggPathway('<?php echo($projectID); ?>', '<?php echo($pathwayID); ?>');
</script>

</body>
</html>