<?php 
session_start();
include('../includes/dbConnect.php');
include('../amfphp/services/Util.php');
include('../includes/restrictAccess.php');

$mode = "single";
$program = "blastn";
if(isset($_GET['view'])){
 	$view = $_GET['view'];
}
else{
	if(isset($_SESSION['view'])){
		$view = $_SESSION['view'];
	}
	else{
		$view = 'blast';
	}
}

(isset($_GET['from'])) ? $previous = $_GET['from'] : $previous = 'none';
if(isset($_GET['db'])){
	$db = $_GET['db'];
	$query = "SELECT default_type FROM project_directory WHERE projectID ='" . $db . "'";
	$result = mysql_query($query);
	$row = mysql_fetch_array($result);
	if($row['default_type'] == 'NT'){
		$defProg = "blastn";
	}
	else{
		$defProg = "blastp";
	}
}
//see if there is an active job in the session table
(isset($_SESSION['active_job'])) ? $active = $_SESSION['active_job'] : $active = 0;
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php $msg = ($mode == "single") ? "BLAST Against SlimeBase Databases" : " Comparative BLAST Against SlimeBase Databases";
echo($msg); ?></title>
<link href="../css/slimebase_basic.css" rel="stylesheet" type="text/css" />
<link href="../css/blast.css" rel="stylesheet" type="text/css" />
<script language='javascript' src='javascript/seqUtils.js' type='text/javascript'></script>
<script language="javascript" src="javascript/blastFunctions.js" type='text/javascript'></script>
<script language='javascript' type='text/javascript'>
var blastParsed = false;
var view = '<?php echo($view); ?>';
var activeResult = <?php echo($active); ?>;
</script>
<link rel="icon" type="image/png" href="../favicon.png" />
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

        <div id='contentContainer'>
        	<div id='blastSubNav'>
            		<div id='blastNavContainer'>
                    <!--<div class='blastNavLabel'>BLAST:</div>-->
                    <div class='blastNavElement'><a href='#' onclick='setView("blast", view); return(false);'>new</a></div><div class='blastNavElement'><a href='#' onclick='setView("queue", view); '>queue</a></div><div class='blastNavElement'><a href='#' onclick='setView("result", view); '>output</a></div>
                	</div>
             </div>
        	<form name='blastOptionsForm' enctype="multipart/form-data"  id='blastOptionsForm' method='post' action='blast/prepareBLAST.php'>
            <div id='leftContent'>
            <div id='blastOptions'>
        		<div id='blastTop'></div>
                <div id='blastMiddle'>
                	<div id='blastOptionContainer'>
                    	<div class='blastOptionElement'><h1>BLAST Program</h1>
                    	<select name='program' id='programSelect' onchange="reloadDatabases(this.value, '<?php echo($mode); ?>');">
                           	<?php if(isset($defProg)){
								if($defProg == 'blastn'){ ?>
                             <option value='blastn' selected="selected">blastn (DNA vs DNA)</option>
                             <option value='blastp'>blastp (protein vs protein)</option>
                             	<?php }
								else{ ?>
                             <option value='blastn'>blastn(DNA vs DNA)</option>
                             <option value='blastp' selected="selected">blastp (protein vs protein)</option>
                             	<?php }}else{
							?>
                      
                            <option value='blastn'>blastn (DNA vs DNA)</option>
                        	<option value='blastp'>blastp (protein vs protein)</option>
                            <?php } ?>
                            <option value='blastx'>blastx (DNA vs protein)</option>
                        	<option value='tblastn'>tblastn (protein vs transl. DNA)</option>
                        	<option value='tblastx'>tblastx (transl. DNA vs transl DNA)</option>
                        </select>
                        </div>
                        <div class='blastOptionElement'>
                        	<h1>Mode</h1>
                            <select name='mode' id='modeSelect'>
                            	<option value='normal'>normal</option>
                                <option value='small'>very short query</option>
                            </select>
                        </div>
                        <div class='blastOptionElement' id='databaseSelectDiv'>
                        	<?php if($mode == 'single'){ ?>
                            <h1>Database</h1>
                            <select name='database' id='databaseSelect'>
                            	<option value='-1'>Select a Database</option>
                            	<?php
									$toEcho = "";
									if($program == "blastn" || $program == "tblastn" || $program == "tblastx"){
										if($universalAccess == false){
											if($groupMembership != ""){
												$query = "SELECT projectID, project_name FROM project_directory, project_group_membership, user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "' AND user_group_membership.group_id = project_group_membership.group_id AND project_directory.projectID = project_group_membership.project_id";
											}
											else{
												$query = "SELECT projectID, project_name FROM project_directory, users, permissions WHERE num_NT_seqs != '0' AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
											}
										}
										else{
											$query = "SELECT projectID, project_name FROM project_directory WHERE num_NT_seqs != '0' ORDER BY project_name";
										}
									}
									else{
										if($universalAccess == false){
											$query = "SELECT projectID, project_name FROM project_directory, users, permissions WHERE num_AA_seqs != '0' AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
										}
										else{
											$query = "SELECT projectID, project_name FROM project_directory WHERE num_AA_seqs != '0' ORDER BY project_name";
										}
									}
									$result = mysql_query($query);
									while($row = mysql_fetch_array($result)){
										if(isset($db)){
											if($db == $row['projectID']){
												$toEcho .= "<option value='" . $row['projectID'] . "' selected>" . $row['project_name'] . "</option>"; 
											}
											else{
												$toEcho = $toEcho . "<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>"; 
											}
										}
										else{
											$toEcho = $toEcho . "<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>"; 
										}
									}
									echo($toEcho);
								?>
                            </select>
							<?php } else{ ?>
                            <h1>Databases</h1>
                            <?php } ?>
                        </div>
                        <div class='blastOptionElement'>
                        	<h1>E-value</h1>
                            <select name='eval' id='evalSelect'>
                            	<option value='1000'>1000</option>
                                <option value='100'>100</option>
                                <option value='10' selected="selected">10</option>
                                <option value='1'>1</option>
                                <option value='1e-04'>1e-04</option>
                                <option value='1e-10'>1e-10</option>
                                <option value='1e-20'>1e-20</option>    
                                <option value='1e-40'>1e-40</option>
                                <option value='1e-50'>1e-50</option>
                            </select>
                        </div>
                        <div class='blastOptionElement'>
                        	<h1>Email (optional)</h1>
                            <input type='text' id='email' name='email' onfocus='this.select();'/>
                        </div>
                       
                    </div>
           	 	</div>
                <div id='blastBottom'></div>	
        	</div>
             <div id='blastAdditional'></div>
            </div>
            <div id='rightContent'>
            <div id='blastMain'>
            	<div id='blastMainTopLeft'></div>
                <div id='blastMainTopMiddle'></div>
                <div id='blastMainTopRight'></div>
                <div id='blastMainMiddle'>
                	<div id='blastSubmit'>
                	<div id='queryInsert'>
                    	<div class='blastHeader'><h1>New BLAST Job</h1></div>
                        <div id='queryDiv'>
                    	<h2>paste sequences</h2>
                        <textarea id='queryEntry' name='queryEntry'></textarea><br /><br />
                        <h2>or upload file</h2>
                        <div id='fileUpload'><input type='file' name='queryFile' id='queryFile' /><a href="javascript:;" onclick='clearInputBox();'>clear</a></div>
                        <br />
                        </div>
                    </div>
                    </div>
                    <div id='blastQueue'>
                    	<div id='blastQueueHeader'><h1>BLAST Results</h1></div>
                    	<div class='blastQueueSubHeader'><h1>Your Finished Jobs</h1></div>
                        <div id='finishedContainer'>
                        <?php $query = "SELECT * FROM finished_jobs WHERE session_id = '" . session_id() . "'";
						$result = mysql_query($query);
						if(mysql_num_rows($result)){
							$counter = 1;
							while($row = mysql_fetch_array($result)){
								echo("<div class='blastJob'>");
								echo("<div class='blastJobCount'><b>" . $counter . ".</b></div>"); 
								echo("<div class='blastJobDescription'><a href='#' onclick=\"showResult(" . $row['blast_id'] . ", 'queue');return(false);\">" . $row['description'] . "</a></div>");
								echo("<div class='blastJobTime'>" . $row['time_finished'] . "</div>");
								echo("</div>");
								$counter++;
							}
						}
						else{
							echo("<div class='noBlastJobDiv'><b>No finished BLAST jobs.</b></div>");
						}
                    ?>
                    	</div>
                    	<div class='blastQueueSubHeader'><h1>Global Queue</h1></div>
                          <div id='blastQueueContainer'>
						  <?php 
							$query = "SELECT project_directory.project_name, blast_queue.program, blast_queue.job_submit, blast_queue.session_id FROM blast_queue, project_directory WHERE blast_queue.resource_id = project_directory.projectID";
							$result = mysql_query($query);
							if(mysql_num_rows($result)){
								$counter = 1;
								echo("<div class='blastQueueDescription' style='margin-bottom: 20px;'>This is the global SlimeBase queue. Your BLAST jobs will appear in <span style='color:#008B10;'>green</span>. This page is automatically updated every ten seconds.</div>");
								while($row = mysql_fetch_array($result)){
									echo("<div class='blastJob'>");
									if($row['session_id'] == session_id()){
										echo("<div class='blastJobCount' style='color: #008B10;'><b>" . $counter . ".</b></div>");
										echo("<div class='blastJobDescription' style='color: #008B10;'>" . $row['program'] . " against " . $row['project_name'] . "</div>");
										echo("<div class='blastJobTime style='color: #008B10;''>" . $row['job_submit'] . "</div>");
									}
									else{
										echo("<div class='blastJobCount'><b>" . $counter . ".</b></div>");
										echo("<div class='blastJobDescription'>" . $row['program'] . " against " . $row['project_name'] . "</div>");
										echo("<div class='blastJobTime'>" . $row['job_submit'] . "</div>");
									}
									echo("</div>");
									$counter++;
								}
							}
							else{
								echo("<div class='noBlastJobDiv'><b>No BLAST jobs in queue.</b></div>");
							}
						?>
                        <script language='javascript' type='text/javascript'>
						var timerID = window.setInterval("reloadQueue()", 10000);
                        </script>
                        </div>
                    </div>
                    <div id='blastResult'>
                    	<div id='activeResultDiv'>
                    	<div class='blastHeader'><h1>BLAST Reports</h1></div>
                        <div id="reportNotLoaded">Oops! It looks like there isn't a BLAST report loaded yet. Please select one of your results. If you haven't performed a BLAST yet, <a href='blast.php?view=blast'>start a new one</a>.</div>
                        <div id='selectResultDiv'>
                        	<form action=''>
                            	<select name='selectResult' id='selectResult' onchange="showResult(this.options[this.selectedIndex].value, '');">
                            		<?php
										$query = "SELECT blast_id, description FROM finished_jobs WHERE session_id ='" . session_id() . "'";
										$result = mysql_query($query);
										if(mysql_num_rows($result)){
											echo("<option value='-1' selected>Choose a result</option>");
											while($row = mysql_fetch_array($result)){
												echo("<option value='" . $row['blast_id'] . "'>" . $row['description'] . "</option>");
											}
										}
										else{
											echo("<option>No Results</option>");
										}
									?>
                            	</select></form>
                        </div> 
                    	</div>
                    </div>
                </div>
                <div id='blastMainBottomLeft'></div>
                <div id='blastMainBottomMiddle'></div>
                <div id='blastMainBottomRight'></div>
            </div>
             <div id='blastStart'><a href="#"><img src='images/blast/blast_start.png' border="0" onclick='checkParams();' alt='Start BLAST Job'/></a></div>
            </div>
            <input type='hidden' name='fromBlast' value='true' />
            </form>
    	</div>
    </div>   
    </div>
</div>
<script language='javascript' type='text/javascript'>if(view != 'default'){setView(view, '<?php echo($previous); ?>');}</script>
<?php	if(isset($_SESSION['loadedResult']) && $view == 'result'){
		echo("<script language='javascript'>showResult(" . $_SESSION['loadedResult'] . ", '');</script>");
	}
?>
</body>
</html>