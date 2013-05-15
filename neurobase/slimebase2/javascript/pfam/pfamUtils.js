var pfamDetailDisplay;
var pfamRoundedBorder = RUZEE.ShadedBorder.create({ corner:10, shadow:20});

function loadPfamPage(id, projectID){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			document.getElementById(selectedPfamPage).removeAttribute("class");
			document.getElementById(id).setAttribute("class", "selected");
			selectedPfamPage = id;
			document.getElementById("pfam-browsing").innerHTML= request.responseText;
			attachPfamRowHandlers(projectID);
		}
	}
	request.open("GET", "ajax/utils/pfam.php?method=loadPage&args=" + projectID + "," + id, true);
	request.send(null);
}

function removeSearchGlass(){
	searchBox.setAttribute("class", "focused-search");
}

function stopRKey(evt) { 
  var evt = (evt) ? evt : ((event) ? event : null); 
  var node = (evt.target) ? evt.target : ((evt.srcElement) ? evt.srcElement : null); 
  if ((evt.keyCode == 13) && (node.id=="pfam-search-element"))  {return false;} 
} 

document.onkeypress = stopRKey;

function returnsCheckSearchText(projectID){
	return function(){
		if(this.value == ''){
			searchBox.setAttribute("class", "idle-search");	
			loadPfamPage(selectedPfamPage, projectID);
		}
		else{
			var request = createAjaxRequest();
			request.onreadystatechange = function(){
				if(request.readyState == 4){
					if(request.responseText != 'null'){
						document.getElementById("pfam-browsing").innerHTML = request.responseText;
						attachPfamRowHandlers(projectID);
					}
				}
			}
			request.open("GET", "ajax/utils/pfam.php?method=updateSearchRes&args=" + projectID + "," + this.value, true);
			request.send(null);
		}
	}
}

function attachPfamRowHandlers(projectID){
	var elements = document.getElementsByClassName('pfam-table-light');
	doAttach(elements, projectID);
	elements = document.getElementsByClassName('pfam-table-dark');
	doAttach(elements, projectID);
}

function doAttach(elements, projectID){
	for(i = 0; i < elements.length; i++){
		elements[i].onmouseover= function(){
			this.style.cursor='pointer';	
			this.style.color="#971b6a";
		};
		elements[i].onmouseout = function(){
			this.style.color="#000000";	
		};
		elements[i].onclick = function(e){
			openPfamDetail(e, this.id, projectID);
			//createPfamDetail(event, this.id, projectID);
		};
	}
}

function openPfamDetail(e, pfamID, projectID){
	if(!e){
		e = event;	
	}
	var dimensions = determineScreenDimensions();
	var leftPos = dimensions[0] / 2 - document.getElementById('pfam-content').offsetWidth/2;
	var topPos = e.pageY;
	window.open("pfam_detail.php?projectID=" + projectID +"&pfamID=" + pfamID, "pfam-details", "resizable=yes,scrollbars=yes,menubar=yes,height=500,width=1210,top=200,left=" + leftPos);
}

function createPfamDetail(event, id, projectID){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var parent = document.getElementById("pfam-browsing");
			if(!(document.getElementById("pfam-detail-container"))){
				pfamDetailDisplay = document.createElement('div');	
				parent.appendChild(pfamDetailDisplay);	
				pfamDetailDisplay.setAttribute("id", "pfam-detail-container");
			}
			var dimensions = determineScreenDimensions();
			var leftPos = dimensions[0] / 2 - document.getElementById('pfam-content').offsetWidth/2;
			var topPos = event.pageY;
			pfamDetailDisplay.innerHTML = request.responseText;
			pfamRoundedBorder.render('pfam-detail-container');
			pfamDetailDisplay.style.position = 'absolute';
			pfamDetailDisplay.style.display = 'block';
			pfamDetailDisplay.style.left = leftPos + 50 + 'px';
			pfamDetailDisplay.style.top = topPos + 30 + 'px';
			//execute newly added scripts within in the response text <- method is in seqUtils.js
			
		}
	}
	request.open("GET", "ajax/utils/pfam.php?method=getPfamDetails&args=" + projectID + "," + id, true);
	request.send(null);
	
}

function closePfamDetail(){
		pfamDetailDisplay = document.getElementById("pfam-detail-container");
		pfamDetailDisplay.style.display =' none';
}

function determineScreenDimensions(){
	var winW,winH;
	if (document.body && document.body.offsetWidth) {
 	winW = document.body.offsetWidth;
 	winH = document.body.offsetHeight;
	}
	if (document.compatMode=='CSS1Compat' &&
    	document.documentElement &&
    	document.documentElement.offsetWidth ) {
 		winW = document.documentElement.offsetWidth;
 		winH = document.documentElement.offsetHeight;
	}
	if (window.innerWidth && window.innerHeight) {
 		winW = window.innerWidth;
 		winH = window.innerHeight;
	}
	var results = new Array();
	results[0] = winW;
	results[1] = winH;
	return results;	
}