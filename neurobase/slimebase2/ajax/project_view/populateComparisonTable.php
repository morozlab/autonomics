<?php
	session_start();
        include($_SERVER['DOCUMENT_ROOT'] . '/includes/restrictAccess.php');
        include($_SERVER['DOCUMENT_ROOT'] . '/includes/dbConnect.php');
	
	$page = $_GET['page'];
	$project1 = $_GET['project1'];
	$project2 = $_GET['project2'];
	
	$query = "SELECT COUNT(*) as c FROM cross_comparisons WHERE project_id_1 ='" . $project1 . "' AND project_id_2 ='" . $project2 . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	$numPages = ceil($row['c']/100);
	if($page < 1){
		$page = 1;
	}
	else if($page > $numPages){
		$page = $numPages;	
	}
	$numResults = $row['c'];
	//echo($query . " " . $project1 . " " . $project2);
	$pageOffset = ($page - 1) * 100;				
	//get the two project names
	$query = "SELECT project_name, projectID FROM project_directory WHERE projectID='" . $project1 . "'";
	$result = mysql_query($query);
	$row1 = mysql_fetch_array($result);
	$projName1 = $row1['project_name'];
	$query = "SELECT project_name FROM project_directory WHERE projectID ='" . $project2 . "'";
	$result = mysql_query($query);
	$row2 = mysql_fetch_array($result);
	$projName2 = $row2['project_name'];
	//begin the fetching and displaying of stuff
	$query = "SELECT sb_id_1, sb_id_2, text1, text2, cross_comparisons.alignment,
`cross_comparisons`.evalue, cross_comparisons.identity FROM cross_comparisons WHERE project_id_1='" . $project1 . "' AND project_id_2='" . $project2 . "' AND sort_type ='1' ORDER BY ranking LIMIT " . $pageOffset . ",100";
	$result = mysql_query($query);
	$toEcho = "<div id='comparison-navigation'>";
	$toEcho .= "<div id='comparison-navigation-left'>";
	if($page != 1){
		$toEcho .= "<a href='#' onclick=\"fetchCrossComparison('" . ($page - 1)  . "');return(false);\" style='margin-right: 5px;'>&lt; previous</a> ";	
	}
	else{
		$toEcho .= "&nbsp;";	
	}
	if($page != $numPages){
		$toEcho .= "<a href='#' onclick=\"fetchCrossComparison('" . ($page + 1) . "');return(false);\">next &gt;</a> ";	
	}
	else{
		$toEcho .= "&nbsp;";	
	}
	$toEcho .= "</div>";
	$toEcho .= "<div id='comparison-navigation-right'><form name='pageSelectForm' id='pageSelectForm'>";
	$toEcho .= "showing page ";
	$toEcho .= "<select name='pageSelect' id='pageSelect' onchange='loadSelectedPage();'>";
	for($i = 1; $i <= $numPages; $i++){
		if($i == $page){
			$toEcho .= "<option value='" . $i . "' selected='selected'>" . $i . "</option>";
		}
		else{
			$toEcho .= "<option value='" . $i . "'>" . $i . "</option>";	
		}
	}
	$toEcho .= "</select>";
	$toEcho .= " of " . $numPages . " <span style='font-size: .8em; margin-left:5px; color:#666;'>(" . $numResults . " homologous sequences)</span></form></div>";
	//show pagination + navigation
	/*$toEcho .= "<div id='page-display'>";
	$toEcho .= "showing page " . $page . " of " . $numPages;
	$toEcho .= "</div>";
	*/
	$toEcho .= "</div>";
	$toEcho .= "<div id='comparison-results'>";
	$toEcho .= "<div class='comparison-result-table'>";
	$toEcho .= "<div class='c-result-header'><div class='sbID-cell'>sbID</div><div class='annot-cell'>annotation</div><div class='eval-cell'>e-value</div><div class='iden-cell'>% identity</div><div class='sbID-cell'>sbID</div><div class='annot-cell'>annotation</div></div>";
	while($row = mysql_fetch_array($result)){
		/*
		if(strlen($row['text2']) > 40){
			$row['text2'] = substr($row['text2'], 0, 40) . "...";
		}
		*/
		//get the information for each sequence
		$text1 = $row['text1'];
		$text2 = $row['text2'];
		if(is_null($text1)){
			$text1 = "No annotation";
		}
		if(strlen($text1) > 40){
			$text1 = substr($text1, 0, 40) . "...";
		}
		if(is_null($text2)){
			$text2 = "No annotation";
		}
		if(strlen($text2) > 40){
			$text2 = substr($text2, 0, 40) . "...";
		}
		$toEcho .= "<div class='c-result-data'>"; 
		$toEcho .= "<div class='sbID-cell'><a href='#' onclick='loadComparisonAlignment(\"" . $project1 . "\", \"". $project2 ."\", \"" . $row['sb_id_1'] . "\", \"" . $row['sb_id_2'] . "\");return(false);' title='Click to view comparative alignment'>" . $row['sb_id_1'] .  "</a></div>";
		
		$toEcho .= "<div class='annot-cell'>"  . $text1 . "</div>";
		$toEcho .= "<div class='eval-cell'>" . $row['evalue'] . "</div><div class='iden-cell'>" . $row['identity'] . "</div>";
		$toEcho .= "<div class='sbID-cell'><a href='#' onclick='loadComparisonAlignment(\"" . $project1 . "\", \"". $project2 ."\", \"" . $row['sb_id_1'] . "\", \"" . $row['sb_id_2'] . "\");return(false);' title='Click to view comparative alignment'>" . $row['sb_id_2'] .  "</a></div>";
		$toEcho .= "<div class='annot-cell'>"  . $text2 . "</div>";
		$toEcho .= "<div id='" . $project1 . "-" . $project2 . "-" . $row['sb_id_1'] . "-" . $row['sb_id_2'] . "'class='c-alignment-container' ></div>";
		$toEcho .= "</div>";
	}
	$toEcho .= "</div>";
	$toEcho .= "</div>";
	echo($toEcho);
	
?>