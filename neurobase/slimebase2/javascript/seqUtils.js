// JavaScript Document

function createAjaxRequest(){
	var request;
	try{
		// Opera 8.0+, Firefox, Safari
		request = new XMLHttpRequest();
	} catch (e){
		// Internet Explorer Browsers
		try{
			request = new ActiveXObject("Msxml2.XMLHTTP");
		} catch (e) {
			try{
			request = new ActiveXObject("Microsoft.XMLHTTP");
			} catch (e){
				// Something went wrong
				alert("Your browser is out of date and does not support Ajax. Please update your browser to view GO Annotations");
				return false;
			}
		}
	}
	return request;
}

function invokeScript(divid)
{
	var scriptObj = divid.getElementsByTagName("SCRIPT");
	var len = scriptObj.length;
	for(var i=0; i<len; i++)
	{
		var scriptText = scriptObj[i].text;
		var scriptFile = scriptObj[i].src
		var scriptTag = document.createElement("SCRIPT");
		if ((scriptFile != null) && (scriptFile != "")){
			scriptTag.src = scriptFile;
		}
		scriptTag.text = scriptText;
		if (!document.getElementsByTagName("HEAD")[0]) {
			document.createElement("HEAD").appendChild(scriptTag)
		}
		else {
			document.getElementsByTagName("HEAD")[0].appendChild(scriptTag);
		}
	}
}

function translateSequence(sequence){
	var request = createAjaxRequest();
	
	request.onreadystatechange = function(){
		if(request.readyState == 4){
		    var resultArray = JSON.parse(request.responseText);
			var toDisplay = "";
			var frame;
			for(var i =0; i < resultArray.length; i++){
				if(i < 3){
					frame = i + 1;
				}
				else{
					frame = (i - 2) * -1;
				}
				toDisplay += "<div class='meat'><h2>Frame " + frame + "</h2></div><div class='translatedSeq'>" + resultArray[i] + "</div>";
			}
			var element = document.getElementById("amino_acid_display");
			element.innerHTML = toDisplay;
			
		}
		
	}
	request.open("GET", "ajax/translateSeq.php?sequence=" + sequence, true);
	request.send(null);
}
