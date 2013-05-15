<?php
session_start();
include_once('../../includes/dbConnect.php');
include('../../includes/restrictAccess.php');
include('../../includes/slimebase2/security/checkPermissions.php');
mysql_select_db('moroz_lab');

$method = $_GET['method'];
$args = $_GET['args'];
$argArray = explode(",", $args);
//first argument is always the projectID
$projectID = $argArray[0];

if($method == "exportPfam"){
	$pfamID = $argArray[1];
	$seqType = $argArray[2];
	$seqSelector;
	if($seqType == "NT"){
		$seqSelector  = "NT";
	}
	else{
		$seqSelector = "AA";
	}
	$toEcho = "";
	$query = "SELECT s1.sb_id, " . $seqSelector . "_sequence, annot, eval, pfamA_acc, evalue, abundance FROM " . $projectID . "_sequences as s1, best_annotations, pfam_annotations WHERE s1.sb_id = pfam_annotations.sb_id AND pfam_annotations.pfamA_acc='" . $pfamID . "' AND best_annotations.project_id=" . $projectID . " AND best_annotations.sb_id = s1.sb_id";
	$result = mysql_query($query);
	$toEcho = "sb_id\tsequence\tbest_annot\tbest_annot_eval\tpfamA_acc\tevalue\tabundance\n";       
	while($row = mysql_fetch_array($result)){
		$toEcho .= $row['sb_id'] . "\t" . $row[$seqSelector . '_sequence'] . "\t" . $row['annot'] . "\t" . $row['eval'] . "\t" . $row['pfamA_acc'] . "\t" . $row['evalue'] . "\t" . $row['abundance'] . "\n";
	}
	header('Content-type: text/plain');
	header("Content-Disposition: attachment; filename=NeuroBase_". $projectID . "_pfam" . $pfamID. "_export.txt");
	echo($toEcho);
}
elseif($method == "exportKEGG"){
	$pathwayID = $argArray[1];
	$seqType = $argArray[2];
	$seqSelector;
	if($seqType == "NT"){
		$seqSelector  = "NT";
	}
	else{
		$seqSelector = "AA";
	}
	$query = "SELECT s1. " . $seqSelector . "_sequence, s1.sb_id, s1.abundance, best_annotations.annot as best_annot, best_annotations.eval as best_annot_eval, kegg_long_id, kegg_pathway_description, kegg_annotations.evalue from " . $projectID . "_sequences as s1, kegg_annotations, kegg_confs, best_annotations WHERE best_annotations.project_id=" . $projectID . " AND best_annotations.sb_id = s1.sb_id AND kegg_annotations.project_id='" . $projectID . "' AND kegg_annotations.path_id='" . $pathwayID . "' AND kegg_annotations.kegg_id = kegg_confs.kegg_id AND kegg_annotations.path_id = kegg_confs.path_id AND kegg_annotations.sb_id = s1.sb_id ORDER BY kegg_annotations.kegg_id";
	$result = mysql_query($query);
	$toEcho = "sb_id\tsequence\tbest_annot\tbest_annot_eval\tkegg_long_id\tkegg_pathway_description\tkegg_evalue\tabundance\n";
	while($row = mysql_fetch_array($result)){
		$toEcho .= $row['sb_id'] . "\t" . $row[$seqSelector . "_sequence"] . "\t" . $row['best_annot'] . "\t" . $row['best_annot_eval'] . "\t" . $row['kegg_long_id'] . "\t" . $row['kegg_pathway_description'] . "\t" . $row['evalue'] . "\t" . $row['abundance'] . "\n";
	}
	header('Content-type: text/plain');
	header("Content-Disposition: attachment; filename=NeuroBase_". $projectID . "_KEGG" . $pathwayID. "_export.txt");
	echo($toEcho);
	
}

?>