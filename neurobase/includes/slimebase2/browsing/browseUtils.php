<?php
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
	}
	else{
		if($projectID != -1){
			//calculate pagination data
			$curPage = 1;
			//check if this is a search or browse
			if(isset($search)){
				//if this does not have the new type of annotation, do this
				if($haveAnnotation == false){
					$query = "SELECT sb_id FROM " . $projectID . "_sequences WHERE description LIKE '%" . $search . "%'";
				}
				else{
					$query = "(SELECT DISTINCT sb_id FROM homology, annotation_db WHERE homology.annotation_id = annotation_db.annotation_id AND annotation_db.text LIKE '%" . $search . "%' AND project_id='" . $projectID . "') UNION (SELECT DISTINCT sb_id FROM " . $projectID . "_sequences WHERE sb_id='" . $search . "')";
				}
				$result = mysql_query($query);
				$numSeqs = mysql_num_rows($result);
			}
			else{
				if($filter == 'none'){
					$query = "SELECT num_AA_seqs, num_NT_seqs, project_name FROM project_directory WHERE projectID ='" . $projectID . "'";
					$result = mysql_query($query);
					$row = mysql_fetch_array($result);
					if($row['num_NT_seqs'] > $row['num_AA_seqs']){
						$numSeqs = $row['num_NT_seqs'];
					}
					else{
						$numSeqs = $row['num_AA_seqs'];
					}
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
					$query = "SELECT COUNT(*) as c FROM (SELECT " . $projectID . "_sequences.sb_id FROM " . $projectID . "_sequences, homology WHERE " . $projectID . "_sequences.sb_id = homology.sb_id AND evalue < '" . $filterValue . "' GROUP BY " . $projectID . "_sequences.sb_id) as t";
					$result = mysql_query($query);
					$row = mysql_fetch_array($result);
					$numSeqs = $row['c'];
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
 
 /***************************************************DISPLAY THE SEARCH INTERFACE****************************************/

function displaySearchOptions($toSearch){
	global $projectID;
	$projectID = $toSearch;
	echo("<div id='resultHeader'>SlimeBase Search</div><div id='pageInfoContainer'></div>");
	displayAddlSearch();
}

function displayProjectOptions(){
	global $projectID;
	global $projectName;
	global $haveGO;
	echo("<div id='optionsWrapper'>");
	echo("<div id='optionLeft'></div>");
	echo("<div id='optionContent'>");
	echo("<div class='optionTitleElement'><a href='http://74.252.103.104:8888/slimebase2/viewSeqs.php?projectID=" . $projectID . "&filter=none' TITLE='Return to first page of project'>" . $projectName . "</a></div>");
	echo("<div id='optionElementsContainer'>");
	echo("<div class='optionElement'>sort</div>");
	//filter menu
	echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>filter</a>");
	echo("<ul><li class='invisible'><br/></li><li class='colored'>by annotation<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=withAnnot';\"><noscript><a href='viewSeqs.php?filter=withAnnot'></noscript>with annotation<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=noAnnot';\"><noscript><a href='viewSeqs.php?filter=noAnnot'></noscript>without annotation<noscript></a></noscript></li></ul></li>");
	echo("<li class='colored'>by evalue<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-04';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-04'></noscript>1e-04<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-10';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-10'></noscript>1e-10<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-20';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-20'></noscript>1e-20<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-40';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-40'></noscript>1e-40<noscript></a></noscript></li></ul></li>");
	echo("<li class='colored'><a href='viewSeqs.php?filter=none'>no filter</a></li>");
	echo("</ul></li></ul></div>");
	//end filter menu ***************************************
	echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>actions</a>");
	echo("<ul><li class='invisible'><br/></li><li class='colored'><a href='viewSeqs.php?searchProj=" . $projectID . "'>search</a></li><li class='colored'><a href='blast.php?db=" . $projectID . "'>BLAST</a></li>");
	if($haveGO == true){
		echo("<li class='colored'><a href='go/GOBrowser.php?projectID=" . $projectID . "'>browse GO</a></li>");
	}
	echo("</ul></li></ul>");
	echo("</div>");
	echo("</div></div>");
	echo("<div id='optionRight'></div>");
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
	//echo("<div id='resultHeader'>" . $projectName . "</div>");
	$pageOffset = ($curPage - 1) * $pageSize;
	$limiter = "LIMIT " . $pageOffset . ", " . $pageSize;
	$orderer = "ORDER BY evalue";
	$seqCounter = $pageOffset + 1;
	//this is where we get the list of sequences from the DB to show
	if($haveAnnotation == true){
		$orderer = "ORDER BY  IFNULL(homology.evalue, 10000)";
		if($filter == "none"){
			$query = "SELECT " . $projectID . "_sequences.sb_id FROM " . $projectID . "_sequences LEFT OUTER JOIN homology ON " . $projectID . "_sequences.sb_id = homology.sb_id GROUP BY " . $projectID . "_sequences.sb_id " . $orderer . " "   . $limiter;
		}
		else if ($filter == "noAnnot"){
			$orderer = "ORDER BY sb_id";
			$query = "SELECT sb_id FROM " . $projectID . "_sequences WHERE NOT EXISTS (SELECT sb_id FROM homology WHERE homology.sb_id = " . $projectID . "_sequences.sb_id AND project_id ='" . $projectID . "') " . $orderer . " " . $limiter;
		}
		else if ($filter == 'withAnnot'){
			$query = "SELECT " . $projectID . "_sequences.sb_id FROM " . $projectID . "_sequences, homology WHERE " . $projectID . "_sequences.sb_id = homology.sb_id GROUP BY " . $projectID . "_sequences.sb_id " . $orderer . " " . $limiter;
		}
		else if ($filter == 'evalue'){
			$query = "SELECT " . $projectID . "_sequences.sb_id FROM " . $projectID . "_sequences, homology WHERE " . $projectID . "_sequences.sb_id = homology.sb_id AND evalue < '" . $filterValue . "' GROUP BY " . $projectID . "_sequences.sb_id " . $orderer . " " . $limiter;
		}
	}
	else{
		$query = "SELECT sb_id, description, length, date FROM " . $projectID . "_sequences " . $limiter;
	}
	$result = mysql_query($query);
	//displays numbers of items currently on page
	echo("<div id='pageInfoContainer'>");
	echo("<div id='resultDescription'>");
	echo("Sequences <b>" . ($pageOffset  + 1) . "</b> - <b>" . ($pageOffset + $pageSize) . "</b> of " . $numSeqs);
	echo("</div>");
	//display the item navigation options
	echo("<div id='resultNav'>");
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
function displaySequence($sb){
	global $display;
	global $projectID;
	//for each sequence in the proejct on this page, do the following
	//outer container
	echo("<div class='sequenceContainer'>");
	//inner container
	echo("<div class='sequenceContainerInner'>");
	if($display == "annotation"){
		//seq header div
		echo("<div class='seqHeaderContainer'>");
		//display the original sequence name as a link
		echo("<div class='sequenceDetailContainer'>");
		echo("<table class='annotationTable' width='100%'>");
		echo("<tr><td class='sbCell'><a href='seqDetail.php?sbid=" . $sb . "&projectID=" . $projectID . "' TITLE='View Sequence Details'>sb|" . $sb . "|</a></td>");
		//display the annotation for this sequence
		$query2 = "SELECT homology.annotation_id, annotation_db.text, homology.evalue, annotation_sources.name FROM homology, annotation_db, annotation_sources WHERE homology.annotation_id = annotation_db.annotation_id AND homology.sb_id = '" . $sb. "' AND annotation_db.source = annotation_sources.id LIMIT 0,1";
		$result2 = mysql_query($query2);
		$firstAnnot = true;
		echo("<td><table>");
		while($row2 = mysql_fetch_array($result2)){
			echo("<tr>");
			$name = $row2['name'];
			if($firstAnnot == true){
				$firstAnnot = false;
			}
			else{
			}
			echo("<td class='iconCell'>");
			//display NR annotation icon
			if($name == "nr"){
			}
			//display SP annotation icon
			else{
				echo("<a href='http://www.uniprot.org/uniprot/" . $row2['annotation_id'] . "' TITLE='View Uniprot Entry'><img border='0' src='images/view_seqs/Swiss-Prot-logo.png' /></a>");
				$externalPath = "http://www.uniprot.org/uniprot/";	
			}
			echo("</td>");
			$description = $row2['text'];
			if(strlen($description) > 77){
				$description = substr($description, 0, 77) . "...";
			}
			echo("<td class='annotationTextCell'>" . $description . "</td>");
			echo("<td class='evalCell'>" . $row2['evalue'] . "</td>");
			echo("</tr>");
		}
		echo("</table></td>");
		if($firstAnnot == true){
			//there we no annotations found, tell the user this sequence has no annotation
			echo("<td>No annotation</td>");
		}
		echo("</tr></table>");
		echo("</div>");
		echo("</div>");
	}
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
	echo("<div id='resultHeader'>Search Results - <span style='color:#000000;'>" . $search . " in <a href='viewSeqs.php?projectID=" . $projectID . "'>" . $projectName . "</a></span></div>");
	echo("<div id='pageInfoContainer'>");
	echo("<div id='resultDescription'>");
	if($numSeqs == 0){
		echo("</div></div><div id='noResults'>No results found for <b>" . $search . "</b>, please try another query.</div>");	
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
			//get list of sequences with annotation matching the search term
			$query = "(SELECT DISTINCT sb_id FROM homology, annotation_db WHERE homology.annotation_id = annotation_db.annotation_id AND annotation_db.text LIKE '%" . $search . "%' AND project_id='" . $projectID . "') UNION (SELECT DISTINCT sb_id FROM " . $projectID . "_sequences WHERE sb_id ='" . $search. "') " . $limiter;
			$result = mysql_query($query);
			//iterate over the found sequences and display them
			while($row = mysql_fetch_array($result)){
				displaySequence($row['sb_id']);
			}
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
	echo("<div id='addlSearch'><div style='float:left;'><form name='addlSearch' id='addlSearch' action='viewSeqs.php' method='get'>");
	echo("search <select name='projectID' id='projectID'>");
	if($universalAccess == false){
					$query = "SELECT projectID, project_name, user_name FROM project_directory, permissions, users WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
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

?>