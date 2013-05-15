<?php
session_start();
$projectID = $_GET['projectID'];
include_once('../../includes/dbConnect.php');
mysql_select_db('moroz_lab');
include('../../includes/restrictAccess.php');
include('../../includes/slimebase2/security/checkPermissions.php');
include('../../includes/slimebase2/resourceCheck.php');
include('../../includes/slimebase2/db/utils.php');


$query = "SELECT project_name, default_type FROM project_directory WHERE projectID='" . $projectID . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$projectName = $row['project_name'];
if($row['default_type'] == "NT"){
	$seqSelector = "NT";	
}
else{
	$seqSelector = "AA";	
}
$query = "SELECT sb_id, description, " . $seqSelector. "_sequence, abundance FROM " . $projectID . "_sequences";
$result = mysql_query($query);
$toEcho = "sb_id\tdescription\tsequence\tabundance\tbest_annotation\tbest_annotation_evalue\tkegg_ids\tkegg_pathways\tpfam_ids\tpfam_descriptions\tgo_ids\tgo_descriptions\n";
while($row = mysql_fetch_array($result)){
	$toEcho .= $row['sb_id'] . "\t" . $row['description'] . "\t" . $row[$seqSelector . '_sequence'] . "\t" . $row['abundance'];
	$query2 = "SELECT annot, eval FROM best_annotations WHERE project_id='" . $projectID . "' AND sb_id=" . $row['sb_id'];
	$result2 = mysql_query($query2);
	if(mysql_num_rows($result2)){
		$annotation = mysql_fetch_array($result2);
		$toEcho .= "\t" . $annotation['annot'] . "\t" . $annotation['eval'];
	}
	else{
		$toEcho .= "\t\t";
	}
	//get the KEGG pathways for this sequence
	if($haveKEGG){
		$s = select();
		$s->table("kegg_annotations");
		$s->columns(array("kegg_id", "kegg_pathway_description"));
		$s->where(_and("project_id='$projectID'", "sb_id='". $row['sb_id'] . "'"));
		$results = $s->execute();
		$kos = '';
		$kegg_descs = "";
		while($row2 = mysql_fetch_array($results)){
			$kos .= $row2['kegg_id'] . ";";
			$kegg_descs .= $row2['kegg_pathway_description'] . ";";
		}
		$toEcho .= "\t" . $kos . "\t" . $kegg_descs;
		
	}
	else{
		$toEcho .= "\t\t";
	}
	//get the pfam domains for this sequence
	if($havePfam){
		$query = "SELECT pfam_annotations.pfamA_acc, description FROM pfam_annotations, pfama WHERE pfam_annotations.pfamA_acc = pfama.pfamA_acc AND pfam_annotations.project_id ='" . $projectID . "' AND sb_id= '" . $row['sb_id'] . "'";
		$result2 = mysql_query($query);
		$pfamAccs = '';
		$pfamDescs = '';
		while($row2 = mysql_fetch_array($result2)){
			$pfamAccs .= $row2['pfamA_acc']	. ";";
			$pfamDescs .= $row2['description'] . ";";
		}
		$toEcho .= "\t" . $pfamAccs . "\t" . $pfamDescs;
	}
	else{
		$toEcho .= "\t\t";	
	}
	//get the GO terms for this sequence
	if($haveGO){
		$s = select();
		$s->columns(array("go_annotation_new.go_id", "description"));
		$s->table(array("go_annotation_new", "go_catalog"));
		$s->where(_and("go_annotation_new.go_id=go_catalog.go_id", "sb_id='" . $row['sb_id'] . "'"));
		$result2 = $s->execute();
		$goIDs = '';
		$goDescs = '';
		while($row2 = mysql_fetch_array($result2)){
			$goIDs .= $row2['go_id'] . ";";
			$goDescs .= $row2['description'] . ";";
		}
		$toEcho .= "\t" . $goIDs . "\t" . $goDescs;
	}
	else{
		$toEcho .= "\t\t";		
	}
	$toEcho .= "\n";
}
header('Content-type: text/plain');
header("Content-Disposition: attachment; filename=". $projectName ."_project_export.txt");
echo($toEcho);


?>