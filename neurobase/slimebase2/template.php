<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<?php 
include('../includes/dbConnect.php');
include('../amfphp/services/Util.php');

?>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title></title>
<link href="../css/slimebase_basic.css" rel="stylesheet" type="text/css" />
</head>

<body>

<div id='container'>
	<div id="content">
	<div id='header'>
    	<img style='float:left; margin-right: 20px;' src="images/slimeicon.jpg" alt="SlimeBase" />
        <div id='navElement'><a href='main.php'>Home</a></div>
        <div id='navElement'><a href='browse.php'>Projects</a></div>
        <div id='navElement'>BLAST</div>
        <div id='navElement'>Quantify</div>
        <div id='navElement'>FASTA</div>
        
    	<div id="search">
        	<form name='search' id='searchForm' method='get' action='viewSeqs.php'>
            	<input type='text' size='20' name='query' value='search...' onFocus='this.select();' />
                <input type='hidden' name='ref' value='browse.php' />
            	<select name='projectID' style='width: 150px;'>
                <option value='-1' selected="selected">Select Database</option>
				<?php
					$query = "SELECT project_name, projectID FROM project_directory WHERE (num_AA_seqs != 0) OR (num_NT_seqs != 0)";
					$result = mysql_query($query);
					while($row = mysql_fetch_array($result)){
						echo("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
					}
				?>
                </select>
            </form>
        </div>
    </div>
    <div id='mainContent'>
    	<div id='contentHeader'></div>
    </div>   
    </div>
</div>
</body>
</html>

<?php 
