<?php

session_start();
include('../../../includes/restrictAccess.php');
include('../../../includes/dbConnect.php');

mysql_select_db('moroz_lab');

$method = $_GET['method'];
$args = $_GET['args'];

if($method == 'updateCart'){
	$argArray = explode("_", $args);
	$project_id = $argArray[1];
	$sb = $argArray[0];
	$query = "INSERT INTO sequence_cart (session_id, sb_id, project_id) VALUES ('" . session_id() . "', '" . $sb . "', '" . $project_id . "')";
	mysql_query($query);
	//get the new number of sequences for this session
	$query = "SELECT COUNT(*) as C FROM sequence_cart WHERE session_id='" . session_id() . "'";
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		echo($row['C']);
	}
}


else if($method == 'displayNavigation'){
	$projectID = $_GET['args'];	
	$toEcho = "<div class='extraContentNavElement' id='sequence-link'><a href='#'>Sequences</a></div>";
	//check if there is GO annotation
	$query = "SELECT  * FROM go_categories WHERE project_id ='" . $projectID . "' LIMIT 0,1";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$toEcho .= "<div class='extraContentNavElement' id='go-link'><a href='go/GOBrowser.php?projectID=" . $projectID . "' target='_blank'>Gene Ontology</a></div>";
		$_SESSION['hasGO'] = true;
	}
	else{
			$_SESSION['hasGO'] = false;	
	}
	$query = "SELECT sb_id FROM pfam_annotations WHERE project_id='" . $projectID . "' LIMIT 0,1";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$toEcho .= "<div class='extraContentNavElement' id='pfam-link'><a href='#'>Pfam Domains</a></div>";
		$_SESSION['hasPfam'] = true;
	}
	else{
		$_SESSION['hasPfam'] = false;	
	}
	
	//check for KEGG pathways
	$query = "SELECT path_id FROM kegg_annotations WHERE project_id='" . $projectID . "' LIMIT 0,1";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$row = mysql_fetch_array($result);
		$toEcho .= "<div class='extraContentNavElement' id='kegg-link'><a href='kegg.php?projectID=" . $projectID . "&pathwayID=" . $row['path_id'] . "' target='_blank'>KEGG</a></div>";	
	}
	$toEcho .= "<div class='extraContentNavElement' id='cross-link'><a href='#'>Homology Search</a></div>";
	$toEcho .= "<div class='extraContentNavElement' id='information-link' ><a href='#'>Information</a></div>";
	echo($toEcho);
}


else if($method == 'loadProjectDetails'){
	$projectID = $_GET['args'];
	$_SESSION['project_loaded'] = $projectID;
	$query = "SELECT project_name, num_NT_seqs, num_AA_seqs, last_mod, num_contigs, num_singlets, average_length, sequencing_tech, browser_img, spec_info_id FROM project_directory WHERE projectID='" . $projectID . "' LIMIT 0,1";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$date = split("-", $row['last_mod']);
	$printDate = $date[1] . "/" . $date[2] . "/" . $date[0];
	//get the number of annotated sequences
	$query = "SELECT COUNT(DISTINCT sb_id) as c FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1'";
	$result2 = mysql_query($query);
	$row2 = mysql_fetch_array($result2);
	$annotated = $row2['c'];
	$numseqs = $row['num_NT_seqs'];
	if($numseqs == 0){
	  $numseqs = $row['num_AA_seqs'];
	}
	$toEcho = $numseqs . " project sequences<br/>";
	//if($row['num_singlets'] != 0){
	//$toEcho .= $row['num_singlets'] . " singlets<br/>";
	//} 
	$toEcho .= $annotated . " annotated sequences<br/>" . "sequenced with: " . $row['sequencing_tech'] . "<br/>last updated on: " . $printDate;
	$echoArray = array();
	$echoArray[0] = $toEcho;
	$echoArray[1] = $row['browser_img'];
	//get the project information for this project
	$query = "SELECT latin_name, common_name, phylum, distribution FROM species_information WHERE info_id ='" . $row['spec_info_id'] . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$echoArray[2] = "<div class='project-information'><h1>" . $row['latin_name'] . "</h1>";
	$echoArray[2] .= "<p><span class='information-label'>common name: </span><span class='information-data'>" . $row['common_name'] . "</span><br/>";
	$echoArray[2] .= "<span class='information-label'>phylum: </span><span class='information-data'>" . $row['phylum'] . "</span><br/>";
	$echoArray[2] .= "<span class='information-label'>distribution: </span><span class='information-data'>" . $row['distribution'] . "</span>";
	echo(json_encode($echoArray));
}

else if($method == 'loadAlign'){
	$argArray = explode(',', $args);
	//select the proejct names
	$query = "SELECT project_name FROM project_directory WHERE projectID ='" . $argArray[0] . "'";
	$result = mysql_query($query);
	$row1 = mysql_fetch_array($result);
	$query = "SELECT project_name FROM project_directory WHERE projectID ='" . $argArray[1] . "'";
	$result = mysql_query($query);
	$row2 = mysql_fetch_array($result);
	$projName1 = $row1['project_name'];
	$projName2 = $row2['project_name'];
	$toEcho = "<div class='comparison-description'>Comparison of:<br/>";
	$toEcho .= "<a href='seqDetail.php?sbid=" . $argArray[2] . "&projectID=" . $argArray[0] . "' title='View sequence details' target='_blank'>" . $argArray[2] . "</a> from " . $projName1 . " and <br/>";
	$toEcho .= "<a href='seqDetail.php?sbid=" . $argArray[3] . "&projectID=" . $argArray[1] . "' title='View sequence details' target='_blank'>" . $argArray[3] . "</a> from " . $projName2 . "<br/>";
	$toEcho .= "</div>";
	$toEcho .= "<br/><br/>";
	$query = "SELECT alignment FROM cross_comparisons WHERE project_id_1 ='" . $argArray[0] . "' AND project_id_2 ='" . $argArray[1] . "' AND sb_id_1 ='" . $argArray[2] . "' AND sb_id_2='" . $argArray[3] . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$toEcho .= $row['alignment'];
	echo($toEcho);
	
}
else if($method == 'reloadSecondSelect'){
	$firstProject = $args;
	$query = "SELECT DISTINCT project_id_2, project_name FROM cross_comparisons_list, project_directory WHERE project_id_1 = '" . $firstProject . "' AND project_id_2 = projectID ORDER BY project_name";
	$result = mysql_query($query);
	$toEcho = "<option value='-1'>Select a second tissue/species</option>";
	while($row = mysql_fetch_array($result)){
		$toEcho .= "<option value='" . $row['project_id_2'] . "'>" . $row['project_name'] . "</option>";
	}
	echo($toEcho);
		
}
else if($method == 'loadExtraContent'){
	$projectID = $args;
	$echoArray = array("", "", "", "");
	$optionList1 = "";
	$optionList2 = "";
	//select all of the projects that have cross-comparison data
	if(isset($_SESSION['group_user'])){
		$query2 = "SELECT DISTINCT project_id_1, project_name FROM cross_comparisons_list, project_directory WHERE project_id_1 = projectID ORDER BY project_name";
	}
	else if($universalAccess == true){
		$query2 = "SELECT DISTINCT project_id_1, project_name FROM cross_comparisons_list, project_directory WHERE project_id_1 = projectID ORDER BY project_name";	
	}
	else{
		$query2 = "SELECT DISTINCT project_id_1, project_name FROM cross_comparisons_list, project_directory, permissions as p1 WHERE project_id_1=projectID AND projectID=p1.project_id AND p1.user_id='" . $_SESSION['user_id'] . "' AND (SELECT user_id FROM permissions WHERE user_id='" . $_SESSION['user_id'] . "' AND project_id=project_id_2)";	
	}
	$result2 = mysql_query($query2);
	$selectedProj = 0;
	$firstProj = true;
	$optionPrefix = "";
	while($row = mysql_fetch_array($result2)){
		if($row['project_id_1'] == $projectID){
			$selectedProj = $row['project_id_1'];
			$optionPrefix = "<option value='" . $row['project_id_1'] . "'>" . $row['project_name'] . "</option>";
		}
		else{
			if($firstProj == true){
				$selectedProj = $row['project_id_1'];
			}
			$optionList1 .= "<option value='" . $row['project_id_1'] . "'>" . $row['project_name'] . "</option>";
		}
		$firstProj = false;
		/*if($row['project_id_1'] != $args){
			//record the first project that isn't the selected project, if we don't have a project selected already
			if($selectedProj == 0){
				$selectedProj = $row['project_id_1'];	
				//$selectedIndex = $index;
				$optionList1 .= "<option value='" . $row['project_id_1'] . "' selected='selected'>" . $row['project_name'] . "</option>";
			}
			else{
				$optionList1 .= "<option value='" . $row['project_id_1'] . "'>" . $row['project_name'] . "</option>";	
			}
			
		}
		else{
			$selectedProj = $row['project_id_1'];
			//$selectedIndex = $index;
			$optionList1 .= "<option value='" . $row['project_id_1'] . "' selected='selected'>" . $row['project_name'] . "</option>";
		}
		//$index++;*/
	}
	$optionList1 = $optionPrefix . $optionList1;
	//select all of the projects that have been compared to the selected 1st project
	$query2 = "SELECT DISTINCT project_id_2, project_name FROM cross_comparisons_list, project_directory WHERE project_id_1 = '" . $selectedProj . "' AND project_id_2 = projectID";
	$result2 = mysql_query($query2);
	while($row = mysql_fetch_array($result2)){
		$optionList2 .= "<option value='" . $row['project_id_2'] . "'>" . $row['project_name'] . "</option>";	
	}
	$echoArray[0] .= "<div id='cross-comparison-selects'><form id='comparison-form' name='comparison-form' name='comparison-form'>compare <select name='comparison-select-1' id='comparison-select-1' class='comparison-select' onchange='reloadSecondSelect();'>";
	$echoArray[0] .= $optionList1;
	$echoArray[0] .= "</select> to <span id='comparison-select-2-span'><select name='comparison-select-2' id='comparison-select-2' class='comparison-select' onchange='fetchCrossComparison(1)'>";
	$echoArray[0] .= $optionList2; 
	$echoArray[0] .= "</select></span></form></div>";		
	$echoArray[1] .= "<div id='information-container'><div class='cross-comparison-header'>This tab contains helpful information such as contact details and organism descriptions.</div><div id='information-content'>";
	$query = "SELECT long_info, latin_name FROM project_information WHERE project_id='" . $projectID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$echoArray[1] .= "<h1>About " . $row['latin_name'] . "</h1>";
	$echoArray[1] .= "<p>" . $row['long_info'] . "</p>";
	//check if there is an abstract
	$query = "SELECT * FROM project_group_membership, abstracts WHERE project_id = '" . $projectID . "' and project_group_membership.group_id = abstracts.group_id LIMIT 0,1";
	$result = mysql_query($query);
	if(mysql_num_rows($result)){
		$row = mysql_fetch_array($result);
		$echoArray[1] .= "<h1>Abstract</h1><span class='abstractTitle' id='abstractTitle'>" . $row['title'] . "</span><p class='abstractAuthors' id='abstractAuthors'>" . $row['authors'] . "</p><p class='abstractText' id='abstractText'>" . $row['text'] . "</p>"; 
	}
	$address = "9505 Oceanshore Blvd<br/>Saint Augustine, FL 32080<br/>";
	$echoArray[1] .= "<h1>NeuroBase was designed and constructed at the University of Florida</h1>";
	$echoArray[1] .= "<div class='people-details'>Professor Leonid L. Moroz<br/>Princple Investigator<br/>" . $address . "904-461-4020<br/><a href='mailto:moroz@whitney.ufl.edu'>moroz@whitney.ufl.edu</a></div>";
	$echoArray[1] .= "<div class='people-details'>Mathew Citarella<br/>Bioinformatics/Design<br/>" . $address . "904-461-4007<br/><a href='mailto:mathew.citarella@gmail.com'>mathew.citarella@gmail.com</a></div>";
	$echoArray[1] .= "<div class='people-details'>Dr. Andrea Kohn<br/>Genomics/Molecular Biology<br/>" . $address . "904-461-4007<br/><a href='mailto:abkohn@msn.com'>abkohn@msn.com</a></div>";
	$echoArray[1] .= "<div class='contact-message'><p>Please direct comments or questions about the database to <a href='mailto:morozgenomics@whitney.ufl.edu'>morozgenomics@whitney.ufl.edu</a>.</p></div>";
	$echoArray[1] .= "</div></div>";
	//load data from pfam
	if($_SESSION['hasPfam'] == true){
		$echoArray[2] = "<div id='pfam-container' class='extraContentContainer'><div class='cross-comparison-header'>This tab contains a summary of all the Pfam domains found in this sequence set.</div><div id='pfam-content'>";
		$echoArray[2] .= "<div id='pfam-heading'>Pfam Domain Browser</div>";
		$echoArray[2] .= "<div id='pfam-search-box'><form><input type='text' id='pfam-search-element' /> </form></div>";
		/*$query = "SELECT COUNT(DISTINCT pfamA_acc) as c FROM pfam_annotations WHERE project_id='" . $projectID . "'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$count = $row['c'];
		$query = "SELECT COUNT(DISTINCT sb_id) as c FROM pfam_annotations WHERE project_id='" . $projectID . "'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$numSeqs = $row['c'];
		$echoArray[2] .= "<div id='pfam-domain-statistics'>" . $count. " domain annotations in " . $numSeqs . " sequences</div>";
		*/
		$echoArray[2] .= "<div id='pfam-navigation'>";
		for ($i=65; $i <= 90; $i++) {
 			$x = chr($i);
			if($x == 'A'){
				$echoArray[2] .= "<a href='#' onclick='loadPfamPage(this.id, $projectID);return(false);' id='$x' class='selected'>$x&nbsp;&nbsp;</a>";
			}
			else{
 				$echoArray[2] .= "<a href='#' onclick='loadPfamPage(this.id, $projectID);return(false);' id='$x'>$x&nbsp;&nbsp;</a>";
			}
		}
		$echoArray[2] .= "</div>";
		$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, pfama.description, pfama.type, counts FROM pfama, pfam_domain_counts WHERE pfama.pfamA_acc = pfam_domain_counts.acc AND project_id='" . $projectID . "' AND pfama.pfamA_id LIKE 'A%'";
		$result = mysql_query($query);
		if(mysql_num_rows($result)){
			$echoArray[2] .= "<div id='pfam-browsing'><table id='pfam-browsing-table'>";
			$echoArray[2] .= "<tr class='pfam-table-header'><td>id</td><td>accession</td><td>description</td><td class='type-column'>type</td><td># of domains</td></tr>";
			$colorator = 0;
			while($row = mysql_fetch_array($result)){
				if($colorator == 0){
					$echoArray[2] .= "<tr class='pfam-table-light' id=\"" . $row['pfamA_acc'] . "\"><td>" . $row['pfamA_id'] . "</td><td>" . $row['pfamA_acc'] . "</td><td>" . $row['description'] . "</td><td>" . $row['type'] . "</td><td align='center'>" . $row['counts'] . "</td></tr>";
					$colorator = 1;
				}
				else{
					$echoArray[2] .= "<tr class='pfam-table-dark' id=\"" . $row['pfamA_acc'] . "\"><td>" . $row['pfamA_id'] . "</td><td>" . $row['pfamA_acc'] . "</td><td>" . $row['description'] . "</td><td>" . $row['type'] . "</td><td align='center'>" . $row['counts'] . "</td></tr>";
						$colorator = 0;
				}
			}
			$echoArray[2] .= "</table>";
			$echoArray[2] .= "</div>";
		}
		$echoArray[2] .= "</div></div>";
	}
	echo json_encode($echoArray);
}

?>