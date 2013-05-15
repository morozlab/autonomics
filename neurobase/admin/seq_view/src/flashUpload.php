<?php
	include("../../../includes/programPaths.php");
	$filename = time() . "_" . $_FILES['flashFile']['name'];
	move_uploaded_file($_FILES['flashFile']['tmp_name'], $tmpDir . "\\" . $filename);	
	echo($filename);
?>