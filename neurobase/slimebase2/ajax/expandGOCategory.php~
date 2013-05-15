<?php
session_start();
include('../../includes/dbConnect.php');
include('../../includes/restrictAccess.php');

$projectID = $_GET['projectID'];
if(isset($_GET['level'])){
	$level = $_GET['level'];
}
$term = $_GET['term'];
if(isset($_GET['expected'])){
	$expected = $_GET['expected'];
}

mysql_select_db('moroz_lab');


$check = "SELECT * FROM go_categories WHERE project_id='" . $projectID . "' LIMIT 0,1";
$checkResult = mysql_query($check);
$returnString = "";
if(!(mysql_num_rows($checkResult))){
if($level == 'top'){
	$check = "SELECT * FROM go_mid_cat WHERE project_id='" . $projectID . "' AND compartment='" . $term . "'";
	$checkResult = mysql_query($check);
	if(!(mysql_num_rows($checkResult))){
	$query = "SELECT go_catalog.go_id, go_catalog.description FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_catalog.go_id = go_annotation.go_higher AND go_catalog.compartment ='" . $term . "' GROUP BY go_catalog.go_id";
	$result = mysql_query($query);
	$total = 0;
	while($row = mysql_fetch_array($result)){
		//for each go_id, count the number we have in the annotation table
		$query = "SELECT COUNT(sb_id) as slim FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_higher ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$slimCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$slimCount = $row2['slim'];
			$total += $slimCount;
		}
		$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"slim\");return(false);'>" . $row['description'] . " (" . $slimCount . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 
	}
	
	if($total != $expected){
	//query the DB for sequences with no go_higher term
	$query = "SELECT DISTINCT(go_catalog.go_id), go_catalog.description FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_catalog.go_id = go_annotation.go_id AND go_catalog.compartment ='" . $term . "' AND go_annotation.go_higher='NA' ORDER BY go_catalog.description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as slim FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$slimCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$slimCount = $row2['slim'];
			$total += $slimCount;
		}
		$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"sequence\");return(false);'>" . $row['description'] . " (" . $slimCount . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 
	}
	}
	
	if($total != $expected){
	//query the DB for seqeunces for which there is no go_higher term in the catalog
	$query = "SELECT DISTINCT(go_catalog.go_id), go_catalog.description FROM go_catalog, go_annotation WHERE go_annotation.project_id ='" . $projectID . "' AND go_catalog.go_id = go_annotation.go_id AND go_catalog.compartment ='" . $term . "' AND go_annotation.go_higher !='NA' AND go_annotation.go_higher NOT IN (SELECT DISTINCT(go_id) FROM go_catalog WHERE compartment='" . $term . "') ORDER BY go_catalog.description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as slim FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$slimCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$slimCount = $row2['slim'];
		}
		$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"sequence\");return(false);'>" . $row['description'] . " (" . $slimCount . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 
	}
	}
	}
	else{
		while($row = mysql_fetch_array($checkResult)){
			$returnString .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"slim\");return(false);'>" . $row['description'] . " (" . $row['count'] . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 		
		}
	}
	
}
else if($level == 'slim'){
	$query = "SELECT DISTINCT(go_catalog.go_id), description FROM go_catalog, go_annotation WHERE go_annotation.go_higher = '" . $term ."' AND go_annotation.go_id = go_catalog.go_id and go_annotation.project_id ='" . $projectID . "' ORDER BY description";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$query = "SELECT COUNT(sb_id) as num FROM go_annotation WHERE project_id = '" . $projectID . "' AND go_id ='" . $row['go_id'] . "'";
		$result2 = mysql_query($query);
		$seqCount = 0;
		while($row2 = mysql_fetch_array($result2)){
			$seqCount = $row2['num'];
		}
		$returnString .= "<div style='float:left; width:840px; padding-bottom:2px; padding-left: 40px;' id='" . $row['go_id'] . "_div'><a href='#' target='_self' id='" . $row['go_id'] . "' onclick='toggleView(this.id, \"expand\", \"sequence\"); return(false);'>" . $row['description'] . " (" . $seqCount . ")</a><div id='" . $row['go_id'] . "_children'></div></div>"; 
	}
}

else if($level == 'sequence'){
	$query = "SELECT go_annotation.sb_id, go_annotation.go_id FROM go_annotation WHERE go_annotation.project_id = '" . $projectID . "' AND go_annotation.go_id ='" . $term . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$sb = $row['sb_id'];
		$query2 = "SELECT text, evalue FROM homology, annotation_db WHERE homology.sb_id = '" . $sb ."' AND homology.annotation_id = annotation_db.annotation_id ORDER BY evalue LIMIT 0,1";
		$result2 = mysql_query($query2);
		$found = false;
		while($row2 = mysql_fetch_array($result2)){
			$found = true;
		$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| " . $row2['text'] . " evalue: " . $row2['evalue'] . "</a></div>"; 
		}
		if($found == false){
			$returnString .= "<div style='float:left; width:740px; padding-bottom:2px; padding-left: 60px;' id='" . $row['go_id'] . "_div'><a style='color:#4969D2;' target='_blank' href='../seqDetail.php?sbid=" . $row['sb_id'] . "&projectID=" . $projectID . "' id='" . $row['go_id'] . "'>" . $row['sb_id'] . "| No Annotation</a></div>";
		}
	}
}
}
else{
	//first, get all the pertinent information for the children of this term
	$query = "SELECT term.acc, term.name, unique_seqs, unique_annotations, term.id FROM term, term2term, go_categories WHERE term2term.term1_id='" . $term . "' AND term.id = term2_id AND go_categories.project_id='" . $projectID . "' AND go_categories.id = term2term.term2_id ORDER BY term.name";
	$result = mysql_query($query);
	$child_categories = false;
	$sub_cats = "";
	while($row = mysql_fetch_array($result)){
		$goID = $row['id'];
		$child_categories = true;
		$sub_cats .= "<div style='float:left; width:860px; padding-bottom:2px; padding-left: 20px;' id='" . $goID . "_div'><a href='#' target='_self' id='" . $goID . "' onclick='toggleCategory(this.id, \"expand\", \"slim\");return(false);'>" . $row['name'] . " (" . $row['unique_seqs'] . " | " . $row['unique_annotations'] . ")</a><div id='" . $goID . "_children'></div></div>"; 	
	}
	
	if(!($child_categories)){
		$goID = $term;
	}
	$query = "SELECT name FROM term WHERE id=" . $term;
	$result = mysql_query($query);
	$name = "";
	while($row = mysql_fetch_array($result)){
		$name = $row['name'];	
	}

	//add a link to show sequences in the category
	$get_cat_seqs = "<div class='go_seq_display'><a href='goCategorySeqs.php?projectID=" . $projectID . "&term=" . $goID . "' target='_BLANK' id='" . $goID. "_seq_link' style='color:#036;'>Get sequences in '" . $name . "'</a></div>";
	$returnString .= $get_cat_seqs . $sub_cats;

}
echo($returnString);


?>