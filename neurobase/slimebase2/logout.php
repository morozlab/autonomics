<?php

session_start();

$ref = $_SERVER['HTTP_REFERER'];

session_destroy();

$rest = "slimebase2/browse.php";

if(isset($_SERVER['HTTP_HOST'])){
	$redirect = "http://" . $_SERVER['HTTP_HOST'] . "/";
}
else{
	$redirect = "http://74.252.103.104:8888/";
}

$redirect .= $rest;

header("Location: " . $redirect );
?>
