<?php
session_start();
$_GET['projectID'] = -5;
include_once('../includes/restrictAccess.php');
include_once('ajax/db.php');
include('../includes/slimebase2/security/checkPermissions.php');



?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<script language='javascript' type='text/javascript' src='js/admin.js'></script>
<script language='javascript' type='text/javascript' src='js/util.js'></script>
<title>Zeroclick Pipeline Administration Panel</title>
<link href="css/admin.css" rel="stylesheet" type="text/css" />
</head>

<body>
<div id='content-container' class='main-content-float'>
	<div id='configure-header'>
    	<div id='configure-nav'>
        	<div class='configure-nav-element'>autorun.configure</div>
            <div class='configure-nav-element'>system.configure</div>
            <div class='configure-nav-element'>project.create</div>
            <div class='configure-nav-element'>queue.status</div>
        </div>
    </div>
	<div id='new-autorun'>
    	<div id='new-autorun-header'>Autorun Configuration</div>
    	<form id='new-project-details' name='new-project-details' action=''>
		<div id='source-selection'>
        	select data source:
            <select id='data-source-select' name='data-source-select'>
            <option selected='SELECTED' value=-1>Please select a data source.</option>
		
		<?php
			//select all of the possible data sources for the system
			$conn = new MySQLi("localhost", "root", "meow12", "zero_click");
			$query = "SELECT * FROM data_sources WHERE active=1";
			$results = $conn->query($query);
			while($row = $results->fetch_assoc()){
				echo("<option value=" . $row['source_id'] . ">" . $row['host'] . " (" . $row['source_type'] . ")</option>");
			}	
		
		?>

        	</select>
            <input type='hidden' name='project-id' id='project-id' value=-1/>
       </div>
       <div id='data-selection'>
            select run: <select id='input-data-select' name='input-data-select'></select>
       </div>
       <div id='new-project-name'>Project name: <input name='input-project-name' id='input-project-name' type='text' /></div>
       <div id='valid-project-name'></div>
       <div id='project-tasks-header'>Select Tasks to Run:</div>
       <div id='project-tasks-checks'>
       </div>
       <div id='configuration-submit'>save configuration</div><div id='proj-submit-progress'></div>
       </form>
	</div>
	<div id='new-task'></div>
	<div id='configure-system'></div>
</div>

</body>
</html>