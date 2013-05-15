<?php 
session_start();
include('../includes/dbConnect.php');
mysql_select_db('moroz_lab');
include('../includes/restrictAccess.php');
include('../includes/slimebase2/layout/pageHeader.php');

?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Quantify Your Sequences in SlimeBase Projects</title>

<link href="../css/quantify.css" rel="stylesheet" type="text/css" />
<link href="../css/page_header.css" rel="stylesheet" type="text/css" />
<link href="../css/slimebase_white.css" rel="stylesheet" type="text/css" />
<script language='javascript' type="text/javascript" src="javascript/quantification/displayScripts.js"></script>
</head>

<body>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<?php include('../includes/slimebaseNav.php'); ?>
        
    <div id="search">
        	<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        </div>
    </div>
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<?php displayPageOptions("quantify", "Quantification"); ?>
        <form id='quantificationForm' name='quantificationForm' enctype="multipart/form-data" method='post' atcion='quant_results.php'>
        <div id='explanation'>SlimeBase's quantification service parses BLAST results to determine the digital expression of your query sequences in a set of raw sequencing reads.</div>
        <div id='step1'>
        	<div class='stepNumber'><img src="images/quantify/1.png" /></div>
            <div class='stepHeader'>provde query sequences</div>
            <div class='stepOptions'>
            	<div id='step1Choices'>
                	<div id='pasteOptionImg'>
                		<a href='#' onclick='showPasteBox();return(false);'><img src="images/quantify/paste_option.png" border='0'/></a>		
                    </div>
                    <div id='uploadOptionImg'>
                    	<a href='#' onclick='showFileBox();return(false);'><img src='images/quantify/upload_option.png' border='0'/></a
                	</div>
                </div>
                <div id='step1Data'>
                	<div id='fileBox'><input type='file' name='queryFile' id='queryFile'/></div>
                    <div id='pasteBox'><textarea id='queryText' name='queryText' onclick='document.quantificationForm.queryText.value="";'>Paste your FASTA-formatted sequences here</textarea></div>
                </div>
            </div>
        </div>
        <div id='step2'>
        	<div class='stepNumber'><img src="images/quantify/2.png" /></div>
        	<div class='stepHeader'>select quantification options</div>
            <div class='stepOptions'>
            	<div id='projectBox'>
                	<select name='projectSelect' class='blueGraySelect' id='projectSelect'>
                    	<option value='-1'>select a data set</option>
                        <?php
							if($universalAccess == false){
								$query = "SELECT projectID, project_name FROM project_directory, permissions, users WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id AND assembly='n' ORDER BY project_name";
							}
							else{
								$query = "SELECT projectID, project_name FROM project_directory WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND assembly='n' ORDER BY project_name";
							}
							$result = mysql_query($query);
							while($row = mysql_fetch_array($result)){
								echo("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
							}
						?>
                    </select>
                </div>
                <div id='evalueSelect'>
                	<select name='evalue' class='blueGraySelect' id='evalue'>
                    	<option valu='-1'>select evalue</option>
                    	<option value='1e-04'>1e-04</option>
                        <option value='1e-10'>1e-10</option>
                        <option value='1e-20'>1e-20</option>
                        <option value='1e-40'>1e-40</option>
                    </select>
                </div>
                <div id='identitySelect'>
                	<select name='identity' id='identity' class='blueGraySelect'>
                    	<option value='-1'>select % identity</option>
                    	<option value='.95'>95%</option>
                        <option value='.90'>90%</option>
                        <option value='.85'>85%</option>
                        <option value='.80'>80%</option>
                        <option value='.70'>70%</option>
                        <option value='.50'>50%</option>
                        <option value='.10'>10%</option>
                    </select>
                </div>               
            </div>
        </div>
        <div id='step3'>
        	<div id='startButton'><img src='images/quantify/start_button.jpg' alt'start quantification'/></div>
        </div>
        </form>
        
    </div>   
    </div>
</div>
</body>
</html>