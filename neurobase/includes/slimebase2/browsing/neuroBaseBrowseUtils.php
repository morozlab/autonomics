<?php
/************************************************* DETERMINE Number of Sequences*****************************************************/
if($projectID != -1){
	//calculate pagination data
	$curPage = 1;
	//check if this is a search or browse
	if($filter == 'none'){
		$query = "SELECT COUNT(DISTINCT sb_id) as c FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$numSeqs = $row['c'];
	}
	else if($filter == 'noAnnot'){
		$query = "SELECT COUNT(sb_id) as c FROM " . $projectID . "_sequences WHERE NOT EXISTS (SELECT sb_id FROM sorted_homology WHERE sorted_homology.sb_id = " . $projectID . "_sequences.sb_id AND project_id ='" . $projectID . "' AND sorted_homology.sort_id='1')";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$numSeqs = $row['c'];
	}
	else if($filter == 'withAnnot'){
		$query = "SELECT COUNT(*) as c FROM (SELECT " . $projectID . "_sequences.sb_id FROM " . $projectID . "_sequences, homology WHERE " . $projectID . "_sequences.sb_id = homology.sb_id GROUP BY homology.sb_id) as t";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$numSeqs = $row['c'];
	}
	else if($filter == 'evalue'){
		$query = "SELECT COUNT(DISTINCT sb_id) as c FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1' AND evalue < '" . $filterValue . "'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$numSeqs = $row['c'];
	}
	else if($filter == 'annotationSource'){
		$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1' AND source='" . $filterValue. "'";
		$result = mysql_query($query);
		$numSeqs = mysql_num_rows($result);
	}
				
				
	//determine the pagination
	$last = ceil($numSeqs/$pageSize);
	if($curPage < 1){
		$curPage = 1;
	}
	else if($curPage > $last){
		$curPage = $last;
	}		
}

/****************************************************DISPLAY Header While Viewing Sequences******************************/

function displaySeqHeader(){
	global $numSeqs;
	global $filter;
	global $filterValue;
	global $projectID;
	$query = "SELECT project_name FROM project_directory WHERE projectID='" . $projectID . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$toEcho = "<div class='project-name-display'>" . $row['project_name'] . "</div>"; 
	$toEcho .= "<div class='project-option-container'>";
	$toEcho .= "<ul class='optionMenu'><li>options";
	$toEcho .= "<ul><li class='invisible'><br/></li><li class='colored'>sort";
	$toEcho .= "<ul><li class='colored'><a href='#' onclick=\"displaySequences(" . $projectID . ",1,'" . $filter . "','" . $filterValue . "', 2);return(false);\">by abundance</a></li><li class='colored'><a href='#' onclick=\"displaySequences(" . $projectID . ",1,'" . $filter . "','" . $filterValue . "', 1);return(false);\">by evalue</a></li>";
	$toEcho .= "</ul>";
	$toEcho .= "</li>";
	$toEcho .= "<li class='colored'><a href='download/downloadProject.php?projectID=" . $projectID . "' id='export_" . $projectID . "' target='_blank'>export</a></li>";
	$toEcho .= "</ul>";
	$toEcho .= "</li>";
	$toEcho .= "</ul>";
	$toEcho .= "</div>";
	echo($toEcho);
}

 /***************************************************DISPLAY THE SEARCH INTERFACE****************************************/

function displaySearchOptions($toSearch){
	global $projectID;
	$projectID = $toSearch;
	echo("<div id='resultHeader'>SlimeBase Search</div><div id='pageInfoContainer'></div>");
	displayAddlSearch();
}

function displayPagination(){
	global $pageSize;
	global $curPage;
	global $numSeqs; 
	$pageOffset = ($curPage - 1) * $pageSize;
	echo("<div class='pageInfoContainer'>");
	echo("<div class='resultDescription'>");
	echo("Sequences <b>" . ($pageOffset  + 1) . "</b> - <b>" . ($pageOffset + $pageSize) . "</b> of " . $numSeqs);
	echo("</div>");
	//display the item navigation options
	echo("<div class='resultNav'>");
	displayPageNav("browse");
	echo("</div>");
	echo("</div>");

}

function displayColumnHeaders(){
	global $haveAbundance;
	global $filter;
	global $filterValue;
	global $projectID;
	echo("<div id='seq-header-container'>");
	//inner container
	echo("<div class='sequenceContainerInner'>");
		echo("<table class='annotationTable' width='100%'>");
		echo("<tr><td class='sbCell' style='text-align:center;'>ID</td>");
		//display the annotation for this sequence
		echo("<td><table>");
		if($filter != 'noAnnot'){
			echo("<tr>");
			echo("<td class='iconCell'>source");
			echo("</td>");
			$abundanceEcho = "";
			$evalClass = "evalCell";
			if($haveAbundance == true){
				$abundanceEcho = "<td class='abundanceCell' title='Sort By Abundance'><a href='#' onclick=\"displaySequences(". $projectID . ", 1, '" . $filter . "','" . $filterValue . "', 2);return(false);\">abundance</a></td>";
				$descriptionLength = 62;
				$evalClass="evalCellSmaller";
			}
			echo("<td class='annotationTextCell'>annotation text</td>");
			echo($abundanceEcho);
			echo("<td class='". $evalClass . "' title='Sort By Evalue'><a href='#' onclick=\"displaySequences(". $projectID . ", 1, '" . $filter . "','" . $filterValue . "', 1);return(false);\">e-value</a></td>");
			echo("</tr>");
		}
		echo("</table></td>");
		if($filter == 'noAnnot'){
			//there we no annotations found, tell the user this sequence has no annotation
			echo("<td>annotation text</td>");
		}
		echo("</tr></table>");
	echo("</div></div>");
}

/***************************************************Display Sequences While Browsing*************************************/
function displayAllSeqs(){
	global $projectID;
	global $pageSize;
	global $pageNum;
	global $last;
	global $numSeqs;
	global $projectName;
	global $numberAnnotations;
	global $display;
	global $filter;
	global $filterValue;
	global $sortVal;
	
	//echo("<div id='resultHeader'>" . $projectName . "</div>");
	$pageOffset = ($pageNum - 1) * $pageSize;
	$limiter = "LIMIT " . $pageOffset . ", " . $pageSize;
	$orderer = "ORDER BY evalue";
	$seqCounter = $pageOffset + 1;
	$filterMod = "";
	//this is where we get the list of sequences from the DB to show
		$orderer = "ORDER BY ranking";
		if($filter == "none"){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='" . $sortVal . "' " . $orderer . " " . $limiter;
		
		}
		else if ($filter == "noAnnot"){
			$query = "SELECT sb_id FROM " . $projectID . "_sequences WHERE NOT EXISTS (SELECT sb_id FROM sorted_homology WHERE sorted_homology.sb_id = " . $projectID . "_sequences.sb_id AND project_id ='" . $projectID . "' and sort_id='1')" . $limiter;
		}
		else if ($filter == 'withAnnot'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id ='" . $projectID . "' AND sort_id='" . $sortVal . "' GROUP BY sb_id " . $orderer . " " . $limiter;
		}
		else if ($filter == 'evalue'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='" . $sortVal . "' AND evalue < '" . $filterValue . "' " . $orderer . " " . $limiter;
		}
		else if($filter == 'annotationSource'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='" . $sortVal . "' AND source ='" . $filterValue . "' " . $orderer . " " . $limiter;
		}
	$result = mysql_query($query);
	//display the item navigation options
	displayPageNav("browse");
	displayColumnHeaders();
	//iterate over the sequences and display the information for each of them
		echo("<div id='seqDisplay'>");
		echo("<div id='seqDisplayBorder'>");
		while($row = mysql_fetch_array($result)){
			displaySequence($row['sb_id']);
			$seqCounter++;
		}
		echo("</div>");
		echo("</div>");
}

/************************************************HANDLES DISPLAYING AN INDIVIDUAL SEQUENCE ENTRY*******************************/
function displaySequence($sb, $search="none"){
	global $projectID;
	global $filter;
	global $filterValue;
	global $haveAbundance;
	global $sortVal;
	//for each sequence in the proejct on this page, do the following
	//outer container
	echo("<div class='sequenceContainer'>");
	//inner container
	echo("<div class='sequenceContainerInner'>");
		echo("<table class='annotationTable' width='100%'>");
		echo("<tr><td class='sbCell'><a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "' TITLE='View Sequence Details' target='_blank'>sb|" . $sb . "|</a></td>");
		//display the annotation for this sequence
		$dbFilter = "";
		$searchConstraint = "";
		if($filter == "annotationSource"){
			$dbFilter = "sorted_homology.source='" . $filterValue . "' AND";
		}
		if($search != "none"){
			if($sbIDSearch == false){
				$searchConstraint = "MATCH (annotation_db.text) AGAINST ('" . $search . "' IN BOOLEAN MODE) AND";
			}
			$dbFilter = "";
		}
		$query2 = "SELECT sorted_homology.annotation_id, annotation_db.text, sorted_homology.evalue, annotation_sources.name, sorted_homology.abundance FROM sorted_homology, annotation_db, annotation_sources WHERE " . $dbFilter . " project_id='" . $projectID . "' AND sort_id='" . $sortVal . "' AND sorted_homology.sb_id = '" . $sb. "' AND  sorted_homology.annotation_id = annotation_db.annotation_id AND " . $searchConstraint . " sorted_homology.source = annotation_sources.id ORDER BY ranking, sorted_homology.evalue LIMIT 0,1";
		$result2 = mysql_query($query2);
		$firstAnnot = true;
		echo("<td><table>");
		while($row2 = mysql_fetch_array($result2)){
			echo("<tr>");
			$name = $row2['name'];
			if($firstAnnot == true){
				$firstAnnot = false;
			}
			echo("<td class='iconCell'>");
			//display NR annotation icon
			if($name == "nr"){
				echo("<a href='http://www.ncbi.nlm.nih.gov/protein/" . $row2['annotation_id'] . "' TITLE='View NR Entry'><img border='0' src='images/view_seqs/white_ncbi.png' /></a>");
			}
			//display SP annotation icon
			else{
				echo("<a href='http://www.uniprot.org/uniprot/" . $row2['annotation_id'] . "' TITLE='View Uniprot Entry'><img border='0' src='images/view_seqs/Swiss-Prot-logo.png' /></a>");	
			}
			echo("</td>");
			$description = $row2['text'];
			$abundanceEcho = "";
			$descriptionLength = 77;
			$evalClass = "evalCell";
			if($haveAbundance == true){
				$abundanceEcho = "<td class='abundanceCell' title='Sequence Abundance'><div class='abundanceImageContainer'><img src='images/view_seqs/abundance_icon.png' /></div><div class='abundanceAmount'>" . $row2['abundance'] . "</div></td>";
				$descriptionLength = 62;
				$evalClass="evalCellSmaller";
			}
			if(strlen($description) > $descriptionLength){
				$description = substr($description, 0, $descriptionLength) . "...";
			}
			echo("<td class='annotationTextCell'>" . $description . "</td>");
			echo($abundanceEcho);
			echo("<td class='". $evalClass . "' title='Annotation E-value'>" . $row2['evalue'] . "</td>");
			echo("</tr>");
		}
		echo("</table></td>");
		if($firstAnnot == true){
			//there we no annotations found, tell the user this sequence has no annotation
			echo("<td>No annotation</td>");
		}
		echo("</tr></table>");
	echo("</div></div>");
}


/************************************************DISPLAY SEQUENCES FROM SEARCH*********************************/

function displaySearchResults(){
	global $projectID;
	global $search;
	global $curPage;
	global $pageSize;
	global $numSeqs;
	global $projectName;
	global $haveAnnotation;
	global $originalSearch;
	global $sort;
	global $sortMod;
	global $sbIDSearch;
	$sortMod = 1;
	if($sort != "none"){
		$query = "SELECT id FROM sort_types WHERE description ='" . $sort . "'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$sortMod = $row['id'];
	}
	else{
		$sortMod = 1;
	}
	echo("<div id='resultHeader'>Search Results - <span style='color:#000000;'>" . $originalSearch . " in <a href='viewSeqs.php?projectID=" . $projectID . "'>" . $projectName . "</a></span></div>");
	echo("<div id='pageInfoContainer1' class='pageInfoContainer'>");
	echo("<div id='resultDescription1' class='resultDescription'>");
	if($numSeqs == 0){
		echo("</div></div><div id='noResults'>No results found for <b>" . $originalSearch . "</b>, please try another query.</div>");	
		displayAddlSearch();	 
	}
	else{
		$pageOffset = ($curPage - 1) * $pageSize;
		//display item numbering
		displayItemNumbering($pageOffset, $pageSize);
		echo("</div>");
		//display the item navigation options
		echo("<div id='resultNav'>");
		displayPageNav("search");
		echo("</div>");
		echo("</div>");
		
		$limiter = "LIMIT " . $pageOffset . ", " . $pageSize;
		if($haveAnnotation == false){
			$limiter = "WHERE description LIKE '%" . $search . "%' " . $limiter;
			$seqCounter = $pageOffset + 1;
			$query = "SELECT sb_id, description, length, date FROM " . $projectID . "_sequences " . $limiter;
			$result = mysql_query($query);
			//iterate over the sequences and display the information for each of them
			echo("<div id='itemEntry'>");
			echo("<table id='seqTable'>");
			while($row = mysql_fetch_array($result)){
				echo("<tr><td id='numberCell'><b>" . $seqCounter . ".</b></td>");
				//split the sb out from the 
				$description = $row['description'];
				$splitDesc = split("\|", $description);
				$sb = $splitDesc[1];
				$description = $splitDesc[2];
				echo("<td id='description'>sb|<a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "'>" . $sb . "</a>|" . $description);
				echo("</td>");
				//reformat the date from the db
				$splitDate = split("\-", $row['date']);
				$date = $splitDate[1] . "/" . $splitDate[2] . "/" . $splitDate[0];
				echo("<td id='seqDate'>");
				echo("</td></tr>");
				$seqCounter++;
			}
			echo("</table>");
			echo("</div>");
		}
		else{
			$seqCounter = $pageOffset + 1;
			if($sbIDSearch == false){
			//get list of sequences with annotation matching the search term
			$query = "SELECT sb_id FROM sorted_homology as h, annotation_db WHERE project_id='" . $projectID . "' AND sort_id='1' AND h.annotation_id = annotation_db.annotation_id AND (MATCH(text) AGAINST ('" . $search . "' IN BOOLEAN MODE)) GROUP BY sb_id ORDER BY ranking " . $limiter;
			}
			else{
				$query = "SELECT sb_id FROM " . $projectID . "_sequences WHERE sb_id='" . $search . "'";
			}
			$result = mysql_query($query);
			//iterate over the found sequences and display them
			echo("<div id='seqDisplay'>");
			echo("<div id='seqDisplayBorder'>");
			while($row = mysql_fetch_array($result)){
				displaySequence($row['sb_id'], $search);
			}
			echo("</div></div>");
		}
	}	
}

/***********************************************DISPLAYS ITEM NUMBERING*********************************************/

function displayItemNumbering($pageOffset, $pageSize){
	global $numSeqs;
	global $last;
	if(($pageOffset + $pageSize) > $numSeqs){
		$endRange = $numSeqs;
	}
	else{
		$endRange = $pageOffset + $pageSize;
	}
	echo("Sequences <b>" . ($pageOffset  + 1) . "</b> - <b>" . $endRange . "</b> of " . $numSeqs);

}


/*************************************************DISPLAYS PAGE NAVIGATION***********************************/
function displayPageNav($mode){
	global $projectID;
	global $pageSize;
	global $last;
	global $pageNum;
	global $numSeqs;
	global $search;
	global $filter;
	global $filterValue;
	global $sortVal;
	//display quick select
	$toEcho = "<div id='sequence-display-header'>&nbsp;";
	$toEcho .= "<div id='quick-view-select'>";
	if($filter == 'none' || $filter == 'withAnnot'){
		$toEcho .= $numSeqs . " annotated sequences | <a href='#' onclick=\"displaySequences(" . $projectID . ", 1, 'noAnnot', 'none', 1);return(false);\">sequences with no annotation</a>";	
	}
	else{
		$toEcho .= "<a href='#' onclick=\"displaySequences(" . $projectID . ",1, 'withAnnot', 'none', 1);return(false);\">sequences with annotation</a> | " . $numSeqs . " unannotated sequences";	
	}
	$toEcho .= "</div>";
	if($pageNum != $last){
		$toEcho .= "<div class='result-link' id='next-link' title='View Page " . ($pageNum + 1) . "'><a href='#' onclick=\"displaySequences(" . $projectID . "," . ($pageNum  + 1). ",'" . $filter . "','" . $filterValue . "', " . $sortVal . ");return(false);\">next ></a></div>";
	}
	//show the page input box
	$toEcho .= "<div id='result-form'><form name='pageSelect' action='viewSeqs.php' method='get'>";
	$toEcho .= "page ";
	$toEcho .= "<select name='curPage' id='seq-page-select' onchange=\"displaySequences(" . $projectID . ", this.options[this.selectedIndex].value, '" . $filter . "', '" . $filterValue . "', " . $sortVal . ");\">";
	for($i = 1; $i <= $last; $i++){
		if($i == $pageNum){
			$toEcho .= "<option value='" . $i . "' selected='selected'>" . $i . "</option>";
		}
		else{
			$toEcho .= "<option value='" . $i . "'>" . $i . "</option>";	
		}
	}
	$toEcho .= "</select>";
	$toEcho .= " of " . $last;
	$toEcho .= "</form></div>"; 
	if($pageNum != 1){
		$toEcho .= "<div class='result-link' title='View Page " . ($pageNum - 1) . "'><a href='#' onclick=\"displaySequences(" . $projectID . "," . ($pageNum - 1). ",'" . $filter . "','" . $filterValue . "', '" . $sortVal . "');return(false);\">< prev</a></div>";	
	}
	$toEcho .= "</div>";
	echo($toEcho);
}

/************************************************DISPLAYS SEARCH AGAIN FORM**************************************/
function displayAddlSearch(){
	global $projectID;
	global $universalAccess;
	global $groupMembership;
	echo("<div id='addlSearch'><div style='float:left;'><form name='addlSearch' id='addlSearch' action='viewSeqs.php' method='get'>");
	echo("search <select name='projectID' id='projectID'>");
	if($universalAccess == false){
		if($groupMembership != ""){
			$query = "SELECT projectID, project_name FROM project_directory, project_group_membership, user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "' AND user_group_membership.group_id = project_group_membership.group_id AND project_directory.projectID = project_group_membership.project_id ORDER BY project_name";
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
		if($row['projectID'] == $projectID){
			echo("<option value='" . $row['projectID'] . "' selected>" . $row['project_name'] . "</option>");
		}
		else{
			echo("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
		}
	}
	echo("</select> for ");
	echo("<input type='text' name='query' /></div><div style='float:left; padding: 5px 0 0 5px;'><a href='#' title='search'><img src='images/search_go_button.png' border=0 onclick='document.addlSearch.submit();' alt='search'/></a></div></form></div>");
	
}


/************************************************DISPLAYS ERROR MESSAGES****************************************/

function displayError($type){
	if($type == 'noDB'){
		echo("<div id='resultHeader'>You must specificy a database to search. Please try again below.</div><div id='pageInfoContainer'></div>");
		displayAddlSearch();
	}
}

/***********************************************FUNCTION TO GET GOANNTOTATION**********************************/

function getGOAnnotation($ID){
	global $projectID;
	$query = "SELECT sb_id FROM go_annotation WHERE sb_id ='" . $ID . "'";
}

function isInteger($n) {

if (preg_match("/[^0-^9]+/",$n) > 0) { return false; } return true;

}

?>