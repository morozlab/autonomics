<?php
include_once('include/std_include.php');
include_once('../config/path.php');

$method = $_GET['method'];
$src = $_GET['src'];

if($method == "getInput"){

	$ret = array();
	$ret['option'] = "";
	//read the contents of the new_runs file for this source
	if(is_numeric($src)){
		$sth = $dbh->prepare("SELECT source_type FROM data_sources WHERE source_id=?");
		$sth->execute(array($src));
		$row = $sth->fetch();
		$ret['type'] = $row['source_type'];
		$ret['option'] = "";
		$sth = $dbh->prepare("SELECT * FROM runname_to_pid WHERE source_id=? AND downloaded='N'");
		$sth->execute(array($src));
		$ret['num_opts'] = $sth->rowCount();
		if($ret['num_opts']){
			while($row = $sth->fetch()){
				$ret['option'] .= "<option value='" . $row['run_name'] . "'>" . $row['run_name'] . "</option>";
			}
		}
		else{
			$ret['option'] = "<option value=-1>No new runs on this data source</option>";	
		}

	}

	echo(json_encode($ret));
	
}



?>