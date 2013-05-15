<?php
include('../../includes/dbConnect.php');

session_start();

mysql_select_db('moroz_lab');

$session = session_id();

$query = "SELECT project_directory.project_name, blast_queue.program, blast_queue.job_submit, blast_queue.session_id FROM blast_queue, project_directory WHERE blast_queue.resource_id = project_directory.projectID";
$result = mysql_query($query);
if(mysql_num_rows($result)){
	$counter = 1;
	echo("<div class='blastQueueDescription' style='margin-bottom: 20px;'>This is the global SlimeBase queue. Your BLAST jobs will appear in <span style='color:#008B10;'>green</span>. This page is automatically updated every ten seconds.</div>");
	while($row = mysql_fetch_array($result)){
		echo("<div class='blastJob'>");
		if($row['session_id'] == $session){
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