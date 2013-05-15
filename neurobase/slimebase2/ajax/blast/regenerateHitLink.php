<?php

include('../../../includes/dbConnect.php');

session_start();

//need to perform security check here

mysql_select_db('moroz_lab');

$changeTo = $_GET['change'];
$resultID = $_GET['resultID'];
$hitID = $_GET['hitID'];
$blastID = $_GET['blastID'];

$query = "SELECT hit_sb FROM active_hits WHERE blast_id ='" . $blastID . "' AND hit_id ='" . $hitID . "' AND result_id ='" . $resultID . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
if($changeTo == 'c'){
	$func = "collapseHit";
}
else if($changeTo == 's'){
	$func = "expandHit";
}
echo("<a href='#' id='" . $resultID . "_" . $hitID . "' onclick=\"" . $func . "(this.id, '" . $blastID . "'); return(false);\" title='Click to expand/collapse hit description'>" . $row['hit_sb'] . "</a>");


?>