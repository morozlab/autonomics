<?php

function recursiveDisplay($parent){
	$query = "SELECT * FROM project_directory WHERE child_of = '" . $parent . "' ORDER BY project_name";
	$result = mysql_query($query);
	if(mysql_num_rows($result) > 0){
		echo('<ul>');
		while($row = mysql_fetch_array($result)){
			$projName = $row['project_name'];
			$projID = $row['projectID'];
			echo('<li>' . $projName);
			recursiveDisplay($projID);
			echo('</li>');
		}
		echo('</ul>');	
	}
}
?>