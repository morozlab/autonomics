<?php


function get_iontorrent_runs($host){
	
	$ret = "";
	$conn = iondb_connect($host);
	echo(pg_errormessage($conn));
	#$results = pg_query($conn, "SELECT fastqLink FROM rundb_results WHERE assembly='N'");
	#while($row = pg_fetch_row($results)){
	#	$fastq = $row['fastqLink'];
#		$name = runname_from_fastqlink($fastq);
#		$ret .= "<option value='$name'>$name</option>";
	#}
	#return $ret;
}

function runname_from_fastqlink($fastq){
	return basename($fastq , ".fastq");		
}



?>