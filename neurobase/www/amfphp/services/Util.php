<?php

class Util{
	
	public function __construct(){
	}
	
	public function translateSequence($params){
		$seqid = $params[0];
		$sequence = $params[1];
		$array = array();
		exec("perl C:\PerlScripts\Format\webTranslate.pl " . $sequence, $array);
		array_push($array, $seqid);
		return $array;	
	}
	
}
?>