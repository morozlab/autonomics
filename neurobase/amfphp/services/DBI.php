<?php
include('../../includes/dbConnect.php');

class DBI{
	
	public function __construct(){
		mysql_select_db("moroz_lab");
	}
	
	public function execSelect($query){
		$ret = array();
		$result = mysql_query($query);
		while($row = mysql_fetch_array($result)){
			array_push($ret, $row);
		}
		return $ret;			
	}
	
	public function execUpdate($query){
		$ret = mysql_query($query);
		return $ret;
	}
}

?>
