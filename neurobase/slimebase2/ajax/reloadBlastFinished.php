<?php
include('../../includes/dbConnect.php');

session_start();

mysql_select_db('moroz_lab');

$session = session_id();

$query = "SELECT * FROM finished_jobs WHERE session_id = '" . $session . "'";
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