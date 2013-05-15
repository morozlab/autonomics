<?php

class FlashDownload{
	
	public function __construct(){
	}
	
	//compreses a project into a zip file and returns the path of that file
	public function compressFile($path){
		include('../../includes/programPaths.php');
		$splitPath = split("\\\\", $path);
		$command = $zipPath . " a -r -tzip " . $splitPath[sizeof($splitPath) - 1] . ".zip " . $path;
		exec($command, $return);
		$oldName =  $wampDocRoot . "amfphp\\services\\" . $splitPath[sizeof($splitPath) - 1] . ".zip";
		$topHalf = "seq_view/tmp/" . $splitPath[sizeof($splitPath) - 1] . time() .  ".zip";
		$newName = $wampDocRoot . $topHalf;
		if(copy($oldName, $newName)){
			unlink($oldName);
			$retPath = "../../../" . $topHalf;
			return	$retPath;
		}
		return "Failure";			
	}
	
	//remove temporary file created during compression
	public function removeTemp($path){
		if(true){
			return 'C:\\wamp\\www\\seq_view\\tmp\\' . $path;
		}
		else{
			return "failure";
		}
	}

}
	
?>