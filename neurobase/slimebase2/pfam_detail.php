<?php
session_start();
include('../includes/dbConnect.php');
include('../includes/restrictAccess.php');
mysql_select_db('moroz_lab');

$projectID = $_GET['projectID'];
$pfamID = $_GET['pfamID'];

include('../includes/slimebase2/security/checkPermissions.php');

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>

<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>NeuroBase Pfam Browser</title>
<link href="../css/pfam.css" rel="stylesheet" type="text/css" />
<link href="../css/pfam_tab.css" rel="stylesheet" type="text/css" />
<link href="../css/prototip.css" rel="stylesheet" type="text/css" />
<!-- load neurobase utility scripts -->
<script type='text/javascript' src="javascript/seqUtils.js"></script>
<script type='text/javascript' src='javascript/displayUtils.js'></script>
<!-- load the prototype library from google --> 
<script type="text/javascript" src="http://www.google.com/jsapi"></script> 
<script type="text/javascript">google.load("prototype", "1.7");</script> 
<!-- load pfam domain libraries -->
<script type='text/javascript' src='javascript/pfam/excanvas.js'></script> 
<script type='text/javascript' src='javascript/pfam/canvas.text.js?reimplement=true&amp;dontUseMoz=true'></script> 
<script type='text/javascript' src='http://pfam.sanger.ac.uk/static/javascripts/faces/optimer-bold-normal.js'></script> 
<script type='text/javascript' src='javascript/pfam/domain_graphics.js'></script> 
<script type='text/javascript' src='javascript/pfam/prettify.js'></script> 
<script type='text/javascript' src="javascript/pfam/prototip.js"></script>
<script type='text/javascript' src="javascript/pfam/styles.js"></script>
<script type='text/javascript' src='javascript/pfam/pfamRendering.js'></script>

</head>

<body>
<div id='pfam-detail-container'>
<?php

	$query = "SELECT pfamA_acc, description, comment FROM pfama WHERE pfamA_acc='" . $pfamID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	//add the pfam description and close button
	$toEcho = "<div id='pfam-detail-header'><a href='http://pfam.sanger.ac.uk/family/" . $pfamID . "' target='_blank'>" . $row['description'] . "(" . $row['pfamA_acc'] . ")</a></div>";
	//add the full HTML info from Pfam (later)
	//add the comment from this entry
	$toEcho .= "<div id='pfam-detail-comment' class='comment'>" . $row['comment'] . "</div>";
	//get the sequences from this project with this pfam domain
	$toEcho .= "<div id='pfam-sequences-header'><span>Sequences with '" . $row['description'] . "' domain </span></div>";
	$toEcho .= "<div id='pfam-sequence-container'>";
	$query = "SELECT " . $projectID . "_sequences.sb_id, annot as best_annot, best_annotations.eval as best_annot_eval, length FROM " . $projectID . "_sequences JOIN pfam_annotations ON " . $projectID . "_sequences.sb_id = pfam_annotations.sb_id LEFT JOIN best_annotations ON best_annotations.sb_id = " . $projectID . "_sequences.sb_id WHERE best_annotations.project_id=" . $projectID . " AND pfam_annotations.project_id=" . $projectID  . " AND pfam_annotations.pfamA_acc='" . $pfamID . "'";
	$result = mysql_query($query);
	$index = 0;
	$toEcho .= "<ol id='pfam-graphics-list'>";
	while($row = mysql_fetch_array($result)){
		//draw the pfam domains for each sequence
		$toEcho .= "<li>";
		$toEcho .= "<div class='pfam-detail-diagram-heading'><a href='seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' target='_blank'>" . $row['sb_id'] . " similar to " . $row['best_annot'] . " (" . $row['best_annot_eval'] . ")</a></div>";
		$toEcho .= "<div class='pfam-graphic-display'><div id='pfam-detail-" . $row['sb_id'] . "-rendering'></div></div>";
		$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, description, evalue, start, stop, type, evalue FROM pfam_annotations, pfama WHERE project_id='" . $projectID . "' AND sb_id='" . $row['sb_id'] . "' AND pfam_annotations.pfamA_acc = pfama.pfamA_acc";	
		$pfamDomainResult = mysql_query($query);
		$toEcho .= "<script language='javascript' type='text/javascript'>featureArray = new Array();";
		while($row2 = mysql_fetch_array($pfamDomainResult)){
			$toEcho .= "featureArray.push(new pfamAnnotation('" . $row2['pfamA_acc'] . "', '" . $row2['start'] . "', '" . $row2['stop'] . "', '" . $row2['description'] . "', '" . $row2['type'] . "', '" . $row2['pfamA_id'] . "', '" . $row2['evalue'] . "', '" . $row2['pfamA_acc'] . "'));\n";	
		}
		$toEcho .= "displayPfamGraphic(featureArray, '" . $row['length'] . "', '" . $row['sb_id'] . "');";
		$toEcho .= "</script></li>";
		$index++;
	}
	$toEcho .= "</ol></div>";
	echo($toEcho);

        //get the sequence type for this project
        $query = "SELECT default_type FROM project_directory WHERE projectID=" . $projectID;
        $result = mysql_query($query);
        while($row = mysql_fetch_array($result)){
          $seq_type = $row['default_type'];
        }
?>
</div>
<div id='pfam-options-container'><?php echo("<div id='pfam-detail-options'><div class='pfam-option'><a href='download/doDownload.php?method=exportPfam&args=" . $projectID . "," . $pfamID . ",". $seq_type . "'>export sequences</a></div></div>"); ?></div>

</body>
</html>