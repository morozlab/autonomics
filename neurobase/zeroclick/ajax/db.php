<?php

$ZC_DSN = "mysql:host=localhost;dbname=zero_click";
$ZC_USER = "root";
$ZC_PASSWD = "meow12";

$dbh = null;
try{
	$dbh = new PDO("mysql:host=localhost;dbname=zero_click", "root", "meow12"); 
	$dbh->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);
}
catch(PDOException $e){
	echo($e->getMessage());
}

function make_db_connection($dsn=null, $user=null, $passwd=null){
	global $ZC_DSN, $ZC_USER, $ZC_PASSWD;
	if(!isset($dns)){
		$dsn = $ZC_DSN;	
	}
	if(!isset($user)){
		$user = $ZC_USER;	
	}
	if(!isset($passwd)){
		$passwd = $ZC_PASSWD;
	}
	return new PDO($dsn, $user, $passwd);
}

function iondb_connect($host){
	
	$conn = null;
	
	$conn = pg_connect("host=10.41.128.76 dbname=iondb user=zeroclick password=Whitney2011");
	echo("moo2");
	return $conn;
		
}


?>