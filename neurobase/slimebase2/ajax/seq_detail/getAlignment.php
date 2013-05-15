<?php
include('../../../includes/dbConnect.php');
mysql_select_db('moroz_lab');
$sb = $_GET['sb'];
$annotationID = $_GET['annotationID'];
$source = $_GET['source'];

$query = "SELECT alignment FROM annotation_alignments WHERE sb_id='" . $sb . "' AND annotation_id ='" . $annotationID . "' AND source='" . $source . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);

echo($row['alignment']);


?>