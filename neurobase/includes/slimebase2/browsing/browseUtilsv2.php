<?php
$resultArray = array();
$resultWithAnnotation = 0;
$idArray = array();
/************************************************* DETERMINE Number of Sequences*****************************************************/
if(!(isset($toSearch))){
	if(isset($_GET['curPage'])){
		$curPage = $_GET['curPage'];
		$numSeqs = $_SESSION['numSeqs'];
		$last = ceil($numSeqs/$pageSize);
		if($curPage < 1){
			$curPage = 1;
		}
		else if($curPage > $last){
			$curPage = $last;
		}
		if(isset($search)){
			getSearchResults();
		}
	}
	else{
		if($projectID != -1){
			//calculate pagination data
			$curPage = 1;
			//check if this is a search or browse
			if(isset($search)){
				getSearchResults();
			}
			else{
				if($haveAnnotation == true){
				if($filter == 'none'){
					$query = "SELECT COUNT(DISTINCT sb_id) as c FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1'";
					$result = mysql_query($query);
					$row = mysql_fetch_array($result);
					$numSeqs = $row['c'];
				}
				else if($filter == 'noAnnot'){
					$query = "SELECT COUNT(sb_id) as c FROM " . $projectID . "_sequences WHERE NOT EXISTS (SELECT sb_id FROM homology WHERE homology.sb_id = " . $projectID . "_sequences.sb_id AND project_id ='" . $projectID . "')";
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
				}
				else{
					$query = "SELECT num_NT_seqs FROM project_directory WHERE projectID='" . $projectID . "'";
					$result = mysql_query($query);
					$row = mysql_fetch_array($result);
					$numSeqs = $row['num_NT_seqs'];
				}
				
			}	
			//determine the pagination
			$last = ceil($numSeqs/$pageSize);
			if($curPage < 1){
				$curPage = 1;
			}
			else if($curPage > $last){
				$curPage = $last;
			}
			$_SESSION['numSeqs'] = $numSeqs;
		}
	}
 }

function getSearchResults(){
	global $originalSearch;
	global $search;
	global $resultArray;
	global $resultWithAnnotation;
	global $numSeqs;
	global $projectID;
	global $haveAnnotation;
	$search = preg_replace("/\s$/", "", $search); 
	$originalSearch = $search;
	$searchArray = explode(" ", $search);
	$newSearch = "";
	$sbIDSearch = false;
	if(sizeof($searchArray) == 1){
		if(isInteger($searchArray[0])){
			$sbIDSearch = true;
		}
		else{
			$newSearch = "+\"" . $searchArray[0] + "\""; //. "*";	
		}
	}
	else{
		for($i = 0; $i < sizeof($searchArray); $i++){
			if(preg_match('/\-/', $searchArray[$i])){
				$searchArray[$i] = "\"" . $searchArray[$i] . "\"";
			}
			if($i == 0){
				$newSearch .= "+" . $searchArray[$i]; 
			}
			else{
				$newSearch .= " +" . $searchArray[$i];
			}
		}
	}
	$resultArray = array();
	if($sbIDSearch == false){
		if($newSearch != ""){
			$search = $newSearch;
		}
		//if this does not have the new type of annotation, do this
		if($haveAnnotation == false){
			$query = "SELECT sb_id, description, date FROM " . $projectID . "_sequences WHERE description LIKE '%" . $search . "%'";
			$result = mysql_query($query);
			$resultArray = generateResultArray($result);
		}
		else{
			//first, get a list of annotation IDs that match the search terms
			$query = "SELECT DISTINCT sb_id FROM (SELECT sorted_homology.annotation_id, sb_id, evalue, text, annotation_db.source, abundance FROM annotation_db, sorted_homology WHERE (MATCH(text) AGAINST ('" . $search . "' IN BOOLEAN MODE)) AND sorted_homology.annotation_id=annotation_db.annotation_id AND project_id='" . $projectID . "' and sort_id='1' ORDER BY ranking) as t1";
			$result = mysql_query($query);
			$resultArray = generateResultArray($result);
			$resultWithAnnotation = 1;			
			
		}
	}
	else{
		$query = "SELECT DISTINCT sb_id FROM " . $projectID . "_sequences WHERE sb_id='" . $search . "'";
		$result = mysql_query($query);
		$resultArray = generateResultArray($result);
	}
	$numSeqs = mysql_num_rows($result);
}
 
function generateResultArray($result){
	$return = array();
	while($row = mysql_fetch_array($result)){
			array_push($return, $row);
	}
	return $return;
}
 
 
 /***************************************************DISPLAY THE SEARCH INTERFACE****************************************/

function displaySearchOptions($toSearch){
	global $projectID;
	$projectID = $toSearch;
	echo("<div id='resultHeader'>SlimeBase Search</div><div id='pageInfoContainer'></div>");
	displayAddlSearch();
}

/***************************************************DISPLAY PROJECT OPTIONS***********************************************/

function displayProjectOptions(){
	global $projectID;
	global $projectName;
	global $haveGO;
	global $haveAbundance;
	global $thisPage;
	echo("<div id='optionsWrapper'>");
	echo("<div id='optionLeft'></div>");
	echo("<div id='optionContent'>");
	if($projectID != -1){
	echo("<div class='optionTitleElement'><a href='http://74.252.103.104:8888/slimebase2/viewSeqs.php?projectID=" . $projectID . "&filter=none' TITLE='Return to first page of project'>" . $projectName . "</a></div>");
	}
	else{
		echo("<div class='optionTitleElement'>" . $projectName . "</div>");
	}
	echo("<div id='optionElementsContainer'>");
	if($haveAbundance == true){
		//sort menu
		echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>sort</a>");
		echo("<ul><li class='invisible'><br/></li><li class='colored'><a href='" . $thisPage . "?sort=evalue'>by evalue</a></li>");
		echo("<li class='colored'><a href='" . $thisPage . "?sort=abundance'>by abundance</a></li>");
		echo("</ul>");
		echo("</div>");
	}
	//filter menu
	echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>filter</a>");
	echo("<ul id='projectOptions'><li class='invisible'><br/></li><li class='colored'>by annotation<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=withAnnot';\"><noscript><a href='viewSeqs.php?filter=none'></noscript>any annotation<noscript></a></noscript></li>");
	$query = "SELECT * FROM annotation_sources"; 
	$result = mysql_query($query);
	while($row = mysql_fetch_array($result)){
		$id = $row['id'];
		$name = $row['name'];
		echo("<li class='colored' onclick=\"location.href='viewSeqs.php?filter=annotationSource&filterValue=" . $id . "';\"><noscript><a href='viewSeqs.php?filter=annotationSource&filterValue=" . $id . "'></noscript>" . $name . " annotation<noscript></a></noscript></li>");
	}
	echo("<li class='colored' onclick=\"location.href='viewSeqs.php?filter=noAnnot';\"><noscript><a href='viewSeqs.php?filter=noAnnot'></noscript>no annotation<noscript></a></noscript></li></ul></li>");
	//start by evalue
	echo("<li class='colored'>by evalue<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-04';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-04'></noscript>1e-04<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-10';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-10'></noscript>1e-10<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-20';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-20'></noscript>1e-20<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-40';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-40'></noscript>1e-40<noscript></a></noscript></li></ul></li>");
	//end by evalue
	echo("<li class='colored'><a href='viewSeqs.php?filter=none'>no filter</a></li>");
	echo("</ul></li></ul></div>");
	//end filter menu ***************************************
	echo("<div class='optionDivider'></div>");
	echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>actions</a>");
	echo("<ul><li class='invisible'><br/></li><li class='colored'><a href='viewSeqs.php?searchProj=" . $projectID . "'>search</a></li><li class='colored'><a href='blast.php?db=" . $projectID . "'>BLAST</a></li>");
	if($haveGO == true){
		echo("<li class='colored'><a href='go/GOBrowser.php?projectID=" . $projectID . "'>browse GO</a></li>");
	}
	echo("</ul></li></ul>");
	echo("</div>");
	echo("</div></div>");
	//echo("<div id='optionRight'>a</div>");
	//echo("<div id='optionFooter'></div>");
	echo("</div>");
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

/***************************************************Display Sequences While Browsing*************************************/
function displayAllSeqs(){
	global $projectID;
	global $pageSize;
	global $curPage;
	global $last;
	global $numSeqs;
	global $projectName;
	global $haveAnnotation;
	global $numberAnnotations;
	global $display;
	global $filter;
	global $filterValue;
	global $sort;
	global $sortMod;
	global $haveAbundance;
	//echo("<div id='resultHeader'>" . $projectName . "</div>");
	$pageOffset = ($curPage - 1) * $pageSize;
	$limiter = "LIMIT " . $pageOffset . ", " . $pageSize;
	$orderer = "ORDER BY evalue";
	$seqCounter = $pageOffset + 1;
	$filterMod = "";
	$sortMod = "";
	if($sort != "none"){
		$query = "SELECT id FROM sort_types WHERE description ='" . $sort . "'";
		$result = mysql_query($query);
		$row = mysql_fetch_array($result);
		$sortMod = $row['id'];
	}
	else{
		$sortMod = 1;
	}
	//this is where we get the list of sequences from the DB to show
	if($haveAnnotation == true){
		$orderer = "ORDER BY ranking";
		if($filter == "none"){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='" . $sortMod . "' " . $orderer . " " . $limiter;
		}
		else if ($filter == "noAnnot"){
			$query = "SELECT sb_id FROM " . $projectID . "_sequences WHERE NOT EXISTS (SELECT sb_id FROM homology WHERE homology.sb_id = " . $projectID . "_sequences.sb_id AND project_id ='" . $projectID . "')" . $limiter;
		}
		else if ($filter == 'withAnnot'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id ='" . $projectID . "' GROUP BY sb_id " . $orderer . " " . $limiter;
		}
		else if ($filter == 'evalue'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1' AND evalue < '" . $filterValue . "' " . $orderer . " " . $limiter;
		}
		else if($filter == 'annotationSource'){
			$query = "SELECT DISTINCT sb_id FROM sorted_homology WHERE project_id='" . $projectID . "' AND sort_id='1' AND source ='" . $filterValue . "' " . $orderer . " " . $limiter;
		}
	}
	else{
		$query = "SELECT sb_id, description, length, date FROM " . $projectID . "_sequences " . $limiter;
	}
	$result = mysql_query($query);
	//displays numbers of items currently on page
	echo("<div id='pageInfoContainer1' class='pageInfoContainer'>");
	echo("<div id='resultDescription1' class='resultDescription'>");
	echo("Sequences <b>" . ($pageOffset  + 1) . "</b> - <b>" . ($pageOffset + $pageSize) . "</b> of " . $numSeqs);
	echo("</div>");
	//display the item navigation options
	echo("<div class='resultNav'>");
	displayPageNav("browse");
	echo("</div>");
	echo("</div>");
	//iterate over the sequences and display the information for each of them
	if($haveAnnotation == true){
		echo("<div id='seqDisplay'>");
		echo("<div id='seqDisplayBorder'>");
		while($row = mysql_fetch_array($result)){
			displaySequence($row['sb_id']);
			$seqCounter++;
		}
		echo("</div></div>");
	}
	else{
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
			//check if the sequence has a GO term
			$annotation = getGOAnnotation($sb);
			$splitDate = split("\-", $row['date']);
			$date = $splitDate[1] . "/" . $splitDate[2] . "/" . $splitDate[0];
			echo("<td id='seqDate'>");
			echo("</td></tr>");
			$seqCounter++;
		}
		echo("</table>");
		echo("</div>");
	}
}

/************************************************HANDLES DISPLAYING AN INDIVIDUAL SEQUENCE ENTRY*******************************/
function displaySequence($sb, $annotationID="none", $search="none"){
	global $display;
	global $projectID;
	global $filter;
	global $filterValue;
	global $haveAbundance;
	global $sort;
	global $sortMod;
	global $sbIDSearch;
	//for each sequence in the proejct on this page, do the following
	//outer container
	echo("<div class='sequenceContainer'>");
	//inner container
	echo("<div class='sequenceContainerInner'>");
	if($display == "annotation"){
		echo("<table class='annotationTable' width='100%'>");
		echo("<tr><td class='sbCell'><a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "' TITLE='View Sequence Details'>sb|" . $sb . "|</a></td>");
		//display the annotation for this sequence
		$dbFilter = "";
		$searchConstraint = "";
		if($filter == "annotationSource"){
			$dbFilter = "sorted_homology.source='" . $filterValue . "' AND";
		}
		if($search != "none"){
			$query2 = "SELECT sorted_homology.annotation_id, annotation_db.text, sorted_homology.evalue, annotation_sources.name, sorted_homology.abundance FROM sorted_homology, annotation_db WHERE project_id='" . $projectID . "' AND sort_id='" . $sortMod . "' AND sorted_homology.sb_id = '" . $sb. "' AND sorted_homology.annnotation_id='" . $annotationID . "' AND sorted_homology.annotation_id = annotation_db.annotation_id ORDER BY ranking LIMIT 0,1";
		}
		else{
			$query2 = "SELECT sorted_homology.annotation_id, annotation_db.text, sorted_homology.evalue, annotation_sources.name, sorted_homology.abundance FROM sorted_homology, annotation_db, annotation_sources WHERE " . $dbFilter . " project_id='" . $projectID . "' AND sort_id='" . $sortMod . "' AND sorted_homology.sb_id = '" . $sb. "' AND  sorted_homology.annotation_id = annotation_db.annotation_id AND " . $searchConstraint . " sorted_homology.source = annotation_sources.id ORDER BY ranking LIMIT 0,1";
		}
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
		}
		echo("</div></div>");
}

function displaySearchResult($result, $resultWithAnnotation=0){
	global $haveAnnotation;
	global $projectID;
	global $originalSearch;
	global $numSeqs;
	global $projectName;
	global $pageOffset;
	global $pageSize;
	global $curPage;
	global $search;
	
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
		$lastIndex = $pageOffset + $pageSize + 1;
		if($lastIndex > sizeof($result)){
			$lastIndex = sizeof($result);
		}
		//case when the project has no annotation
		if($haveAnnotation == false){
			echo("<div id='itemEntry'>");
			echo("<table id='seqTable'>");
			$seqCounter = $pageOffset + 1;
			for($i = $pageOffset; $lastIndex; $i++){
				$description = $result[$i]['description'];
				$date = $result[$i]['date'];
				echo("<tr><td id='numberCell'><b>" . $seqCounter . ".</b></td>");
				//split the sb out from the description
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
		//case when the project has annotation
		else{
			echo("<div id='seqDisplay'>");
			echo("<div id='seqDisplayBorder'>");
			$toEcho = "";
			//don't need to look up annotation information
			/*if($resultWithAnnotation == 1){
				for($i = $pageOffset; $i < $lastIndex; $i++){
					$toEcho .= doDisplaySequence($result[$i]['sb_id'], $result[$i]['annotation_id'], $result[$i]['text'], $result[$i]['evalue'], $result[$i]['source'], $result[$i]['abundance']);	
				}
					
			}
			//need to look up annotation information for each sequence
			else{*/
				for($i = $pageOffset; $i < $lastIndex; $i++){
					$query = "SELECT sorted_homology.annotation_id, text, evalue, abundance, sorted_homology.source FROM sorted_homology, annotation_db WHERE project_id='" . $projectID . "' AND sorted_homology.annotation_id = annotation_db.annotation_id AND sort_id='1' AND sb_id='" . $result[$i]['sb_id'] . "' ORDER BY ranking LIMIT 0,1";
					$result2 = mysql_query($query);
					$row = mysql_fetch_array($result2);
					$toEcho .= doDisplaySequence($result[$i]['sb_id'], $row['annotation_id'], $row['text'], $row['evalue'], $row['source'], $row['abundance']);
				}
			//}
			echo($toEcho);
			echo("</div></div>");
		}
	}
}

function doDisplaySequence($sb, $annotationID='none', $text='none', $evalue='none', $db='none', $abundance='none'){		
	global $display;
	global $projectID;
	global $filter;
	global $filterValue;
	global $haveAbundance;
	global $sbIDSearch;
	$returnString = "";
	//for each sequence in the proejct on this page, do the following
	//outer container
	$returnString .= "<div class='sequenceContainer'>";
	//inner container
	$returnString .= "<div class='sequenceContainerInner'>";
	$returnString .= "<table class='annotationTable' width='100%'>";
	$returnString .= "<tr><td class='sbCell'><a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "' TITLE='View Sequence Details'>sb|" . $sb . "|</a></td>";
		//display the annotation for this sequence
	$firstAnnot = true;
	$returnString .= "<td><table>";
	if($annotationID != "none"){
		$returnString .= "<tr>";
		if($firstAnnot == true){
			$firstAnnot = false;
		}
		$returnString .= "<td class='iconCell'>";
		//display NR annotation icon
		if($db == 1){
			$returnString .= "<a href='http://www.ncbi.nlm.nih.gov/protein/" . $annotationID . "' TITLE='View NR Entry'><img border='0' src='images/view_seqs/white_ncbi.png' /></a>";
		}
		//display SP annotation icon
		else{
			$returnString .= "<a href='http://www.uniprot.org/uniprot/" . $annotationID . "' TITLE='View Uniprot Entry'><img border='0' src='images/view_seqs/Swiss-Prot-logo.png' /></a>";	
		}
		$returnString .= "</td>";
		$abundanceEcho = "";
		$descriptionLength = 77;
		$evalClass = "evalCell";
		if($abundance != "none"){
			$abundanceEcho = "<td class='abundanceCell' title='Sequence Abundance'><div class='abundanceImageContainer'><img src='images/view_seqs/abundance_icon.png' /></div><div class='abundanceAmount'>" . $abundance . "</div></td>";
			$descriptionLength = 62;
			$evalClass="evalCellSmaller";
		}
		if(strlen($text) > $descriptionLength){
			$text = substr($text, 0, $descriptionLength) . "...";
		}
		$returnString .= "<td class='annotationTextCell'>" . $text . "</td>";
		$returnString .= $abundanceEcho;
		$returnString .= "<td class='". $evalClass . "' title='Annotation E-value'>" . $evalue . "</td>";
		$returnString .= "</tr>";
	}
	$returnString .= "</table></td>";
	if($firstAnnot == true){
		//there we no annotations found, tell the user this sequence has no annotation
		$returnString .= "<td>No annotation</td>";
	}
	$returnString .= "</tr></table>";
	$returnString .= "</div></div>";
	return $returnString;
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
	global $idArray;
	global $annotationArray;
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
			//iterate over the found sequences and display them
			echo("<div id='seqDisplay'>");
			echo("<div id='seqDisplayBorder'>");
			foreach($annotationArray as $annotationID => $sbList){
				for($i =0; $i < sizeof($sbList); $i++){
					//displaySequence($sbList[$i], $annotationID, $search);
				}
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
	global $curPage;
	global $numSeqs;
	global $search;
	$forwardString = $projectID . "&pageSize=" . $pageSize . "&curPage=" . ($curPage + 1) . "&numSeqs=" . $numSeqs;
	$reverseString = $projectID . "&pageSize=" . $pageSize . "&curPage=" . ($curPage - 1) . "&numSeqs=" . $numSeqs;
	if($mode == 'search'){
		$forwardString .= "&query=" . $search;
		$reverseString .= "&query=" . $search;
	}
	if($curPage != $last){
		echo("<div id='resultLink'><a href='viewSeqs.php?projectID=" . $forwardString .  "'>next</a></div>");
	}
	//show the page input box
	echo("<div id='resultForm'><form name='pageSelect' action='viewSeqs.php' method='get'>");
	echo("page <input type='text' size='" . strlen($curPage) . "' name='curPage' value='" . $curPage . "' onfocus='this.select();'/> of " . $last);
	echo("<input type='hidden' name='pageSize' value='" . $pageSize . "'/>");
	echo("<input type='hidden' name='projectID' value='" . $projectID . "' />");
	echo("<input type='hidden' name='numSeqs' value='" . $numSeqs . "' />");
	if($mode == 'search'){
		echo("<input type='hidden' name='query' value='" . $search . "' />");
	}
	echo("</form></div>"); 
	//end page input box code
	if($curPage != 1){
		echo("<div id='resultLink'><a href='viewSeqs.php?projectID=" . $reverseString. "'>previous</a></div>");	
	}
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