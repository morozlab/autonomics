<?php


function displayPageOptions($mode, $displayName, $projectID = -1){
	global $haveGO;
	echo("<div id='optionsWrapper'>");
	echo("<div id='optionLeft'></div>");
	echo("<div id='optionContent'>");
	if($mode == 'seqDisplay'){
		if($projectID != -1){
		echo("<div class='optionTitleElement'><a href='http://74.252.103.104:8888/slimebase2/viewSeqs.php?projectID=" . $projectID . "&filter=none' TITLE='Return to first page of project'>" . $displayName . "</a></div>");
		}
		else{
			echo("<div class='optionTitleElement'>" . $projectName . "</div>");
		}
		echo("<div id='optionElementsContainer'>");
		//echo("<div class='optionElement'>sort</div>");
		//filter menu
		echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>filter</a>");
		echo("<ul id='projectOptions'><li class='invisible'><br/></li><li class='colored'>by annotation<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=withAnnot';\"><noscript><a href='viewSeqs.php?filter=withAnnot'></noscript>any annotation<noscript></a></noscript></li>");
		addDBFilters();
		echo("<li class='colored' onclick=\"location.href='viewSeqs.php?filter=noAnnot';\"><noscript><a href='viewSeqs.php?filter=noAnnot'></noscript>no annotation<noscript></a></noscript></li></ul></li>");
	//start by evalue
		echo("<li class='colored'>by evalue<ul><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-04';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-04'></noscript>1e-04<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-10';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-10'></noscript>1e-10<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-20';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-20'></noscript>1e-20<noscript></a></noscript></li><li class='colored' onclick=\"location.href='viewSeqs.php?filter=evalue&filterValue=1e-40';\"><noscript><a href='viewSeqs.php?filter=evalue&filterValue=1e-40'></noscript>1e-40<noscript></a></noscript></li></ul></li>");
	//end by evalue
		echo("<li class='colored'><a href='viewSeqs.php?filter=none'>no filter</a></li>");
		echo("</ul></li></ul></div>");
	//end filter menu ***************************************
		echo("<div class='optionDivider'></div>");
		echo("<div class='optionElement'><ul class='optionMenu'><li><a href='#'>actions</a>");
		echo("<ul><li class='invisible'><br/></li><li class='colored'><a href='viewSeqs.php?searchProj=" . $projectID . "'>search</a></li><li class='colored'><a href='blast.php?db=" . $projectID . "'>BLAST</a></li>");
		if($haveGO == true){
			echo("<li class='colored'><a href='go/GOBrowser.php?projectID=" . $projectID . "'>browse GO</a></li>");
		}
		echo("</ul></li></ul>");
		echo("</div>");
		echo("</div>");
	}
	else if($mode == 'quantify'){
		echo("<div class='optionTitleElement'>" . $displayName . "</div>");
	}
	echo("</div></div>");
	//echo("<div id='optionRight'>a</div>");
	//echo("<div id='optionFooter'></div>");
}

?>