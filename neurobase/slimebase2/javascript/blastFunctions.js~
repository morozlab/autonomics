// JavaScript Document
function reloadDatabases(program, mode){
	if(mode == "single"){
		var el = document.getElementById('databaseSelectDiv');
		var request = createAjaxRequest();
		request.onreadystatechange = function(){
			if(request.readyState == 4){
				var optionText = request.responseText;
				optionText = "<h1>Database</h1><select name='database' id='databaseSelect'><option value='-1'>Select a Database</option>" + optionText + "</select>";
				el.innerHTML = optionText;
			}
		}
		request.open("GET", "ajax/retrieveBlastDBs.php?program=" + program, true);
		request.send(null);
	}
	
}

function checkSession(){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			if(request.responseText == 'no'){
				alert("Your session has expired, redirecting to main blast page");
				window.location.href = 'blast.php';
			}
		}
	}
	request.open("GET", "ajax/checkSession.php", true);
	request.send(null);
}

function reloadQueue(sessionID){
	var el = document.getElementById('blastQueueContainer');
	var request = createAjaxRequest();
	checkSession();
	request.onreadystatechange = function(){
		if(request.readyState  == 4){
			el.innerHTML = request.responseText;
			reloadFinished();
		}
	}
	request.open("GET", "ajax/reloadBlastQueue.php", true);
	request.send(null);
}

function reloadFinished(sessionID){
	var el = document.getElementById('finishedContainer');
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState  == 4){
			el.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/reloadBlastFinished.php", true);
	request.send(null);
}

function checkParams(){
	var el = document.getElementById('blastOptionsForm');
	if(el.database.value == -1){
		alert("Please select a database to BLAST against.");
		return;
	}
	if(el.queryEntry.value == "" && (el.queryFile.value == "")){
		alert("Please paste query sequences or upload a file to BLAST.");
		return;
	}
	document.blastOptionsForm.submit();
}

function clearInputBox(){
	var el = document.getElementById('fileUpload');
	el.innerHTML = el.innerHTML;
}

function jumpResult(index){
	window.location.href = "#result" + index;
}

function setView(current, previous){
	view = current;
	//hide the previous, if there was one
	toggleVisible('blastSubmit', 0);
	toggleVisible('blastStart', 0);
	toggleVisible('blastQueue', 0);
	//reset the main display to the original size	
	makeDisplayDefault();
	toggleVisible('blastResult', 0);
	toggleVisible('leftContent', 0);
	//display the new page section
	if(view == 'queue'){
		toggleVisible('blastQueue', 1);
		toggleVisible('leftContent', 1);
	}
	else if(view == 'default'){
		toggleVisible('leftContent', 1);
		toggleVisible('blastSubmit', 1);
		toggleVisible('blastStart', 1);
	}
	else if(view == 'result'){
		toggleVisible('blastResult', 1);
		//resize the blast display
		var el = document.getElementById('blastMain');
		el.style.width = "900px";
		el.style.marginLeft = "0px";
		el = document.getElementById('blastMainMiddle');
		el.style.width = "868px";
		el = document.getElementById("blastMainTopMiddle");
		el.style.width = "870px";
		el = document.getElementById("blastMainBottomMiddle");
		el.style.width = "870px";
	}
	else if(view == 'blast'){
		toggleVisible('leftContent', 1);
		toggleVisible('blastSubmit', 1);
		toggleVisible('blastStart', 1);
	}
	//set session variablew for view
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
				
		}
	}
	request.open("GET", "ajax/blast/setView.php?view=" + view, true);
	request.send(null);
		
}


function makeDisplayDefault(){
	var el = document.getElementById('blastMain');
	el.style.width = "680px";
	el.style.marginLeft = "20px";
	el = document.getElementById('blastMainMiddle');
	el.style.width = "648px";
	el = document.getElementById("blastMainTopMiddle");
	el.style.width = "650px";
	el = document.getElementById("blastMainBottomMiddle");
	el.style.width = "650px";
}

function toggleVisible(ID, mode){
	var el = document.getElementById(ID);
	if(mode == 1){
		el.style.display = 'block';
	}
	else{
		el.style.display = 'none';
	}
}

function showResult(jobID, from){
	//first, make the BLAST div visible
	setView('result', from);
	if(jobID != activeResult){
		//need to parse and show the appropriate blast result
		var request = createAjaxRequest();
		var el = document.getElementById('blastResult');
		//check if we need to remove a previous node
		while(document.getElementById('activeResultDiv') != null){
			var old = document.getElementById('activeResultDiv');
			old.parentNode.removeChild(old);
		}
		request.onreadystatechange = function(){
			if(request.readyState == 4){
				//only do this if we're still on the blast page
				if(view == 'result'){
					//check if we need to remove a previous node
					while(document.getElementById('activeResultDiv') != null){
						var old = document.getElementById('activeResultDiv');
						old.parentNode.removeChild(old);
					}
					var newDiv = document.createElement("div");
					newDiv.id = "activeResultDiv";
					newDiv.innerHTML = request.responseText;
					el.appendChild(newDiv);
					activeResult = jobID;
				}
			}
		}
		//add a node to show the loading icon
		var div = document.createElement("div");
		div.id = "activeResultDiv";
		div.innerHTML = "<div class='blastHeader'><h1>BLAST Reports</h1></div> <div id='reportNotLoaded'>Loading BLAST report...</div>"
		el.appendChild(div);
		request.open("GET", "ajax/parseBlastResult.php?blastID=" + jobID, true);
		request.send(null);
	}
	
}

function showSequenceDetail(idPart, program, projectID, seqID, blastID){
	var idents = idPart.split("_");
	var resultID = idents[0];
	var hitID = idents[1];
	var container = document.getElementById(idPart + "_full");
	var request = createAjaxRequest();
	container.style.display = 'block';
	container.style.marginLeft = "20px";
	container.style.marginTop = "20px";
	container.style.border = "1px dotted #CCCCCC";
	container.style.padding = "10px";
	container.style.marginBottom = "20px";
	container.innerHTML = "Loading sequence...";
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			container.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/blast/showSequenceDetail.php?seqID=" + seqID + "&program=" + program + "&projectID=" + projectID, true);
	request.send(null);
}

function expandHit(idPart, blastID){
	var idents = idPart.split("_");
	var resultID = idents[0];
	var hitID = idents[1];
	var container = document.getElementById(idPart + "_full");
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			container.style.display = "block";
			container.style.marginLeft = "20px";
			container.style.marginTop = "20px";
			container.style.border = "1px dotted #CCCCCC";
			container.style.padding = "10px";
			container.style.marginBottom = "20px";
			container.innerHTML = request.responseText;
			switchLink("c", resultID, hitID, blastID);
		}
	}
	
	request.open("GET", "ajax/blast/doExpandHit.php?blastID=" + blastID + "&resultID=" + resultID + "&hitID=" + hitID, true);
	request.send(null);
}

function collapseHit(idPart, blastID){
	//all we need to do is hide the hit container, so that when we click again all we need to do is show it instead of
	//parsing/loading everything again
	var container = document.getElementById(idPart + "_full");
	container.style.display = "none";
	var ids = idPart.split("_");
	var resultID = ids[0];
	var hitID = ids[1];
	//change the link
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var linkContainer = document.getElementById(idPart + "_link");
			linkContainer.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/blast/regenerateHitLink.php?change=s&resultID=" + resultID + "&hitID=" + hitID + "&blastID=" + blastID, true);
	request.send(null);
	
}

function showHit(idPart, blastID){
	var container = document.getElementById(idPart + "_full");
	container.style.display = "block";
	var ids = idPart.split("_");
	var resultID = ids[0];
	var hitID = ids[1];
	//change the link
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var linkContainer = document.getElementById(idPart + "_link");
			linkContainer.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/blast/regenerateHitLink.php?change=c&resultID=" + resultID + "&hitID=" + hitID + "&blastID=" + blastID, true);
	request.send(null);
	
}

function switchLink(changeTo, resultID, hitID, blastID){
	//if we need to change the link to collapse, we've already loaded the hit detail
	//this means that the collapse just hides the hit detail, and the next click will show it
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var el = document.getElementById(resultID + "_" + hitID + "_link");
			el.innerHTML = request.responseText;
		}
	}
	
	request.open("GET", "ajax/blast/regenerateHitLink.php?change=" + changeTo + "&resultID=" + resultID + "&hitID=" + hitID + "&blastID=" + blastID, true);
	request.send(null);
}

function checkLength(){
	var inputBox = document.getElementById('queryEntry');
	var smallSeqSelect = document.getElementById('small-sequence-select');
	
}