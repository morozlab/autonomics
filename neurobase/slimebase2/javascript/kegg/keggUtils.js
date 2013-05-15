// JavaScript Document

function drawKeggPathway(projectID, pathwayID){
	//get the image size that is going to be displayed
	var pathImg = new Image();
	var height;
	var width;
	pathImg.onload = function(){
		initCanvas(projectID, pathwayID, pathImg);	
	}
	pathImg.src = 'kegg/map_graphics/map' + pathwayID + '.png';
	
}

function initCanvas(projectID, pathwayID, pathImg){
	//display the base image on the kegg browser
	var canvas = oCanvas.create({
			canvas: "#kegg-graphics",
		});
	canvas.height = pathImg.height;
	canvas.width = pathImg.width;
	
	var displayImg = canvas.display.image({
		x:0,
		y:0,
		image: pathImg.src,
	});
	
	canvas.addChild(displayImg);
	
	//get all of the sequences from this project matching this pathway
	var request = createAjaxRequest();
	
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var returnArray = JSON.parse(request.responseText);
			colorPathComponents(returnArray, canvas);	
			displayPathComponents(returnArray);		
		}
		
	}
	request.open("GET", "ajax/utils/kegg.php?method=colorPath&args=" + projectID + "," + pathwayID, true);
	request.send(null);
	
}

function colorPathComponents(componentArray, canvas){
	var drawn = {};
	for(i = 0; i < componentArray.length; i++){
		var starts = componentArray[i]['start'].split(",");
		var startx = starts[0], starty = starts[1];
		if(startx in drawn){
			continue;
		}
		else{
			drawn[startx] = 1;
		}
		var type = componentArray[i]['type'];
		if(type == "rectangle"){
			var ends = componentArray[i]['end'].split(",");
			var rWidth = ends[0] - startx + 1, rHeight = ends[1] - starty + 1;
			var rectangle = canvas.display.rectangle({
				x: parseInt(startx),
				y: parseInt(starty),
				width: rWidth,
				height: rHeight,
				fill: "#6CDA70",
				opacity: .25
			});
			canvas.addChild(rectangle);
		}
		
		
	}
}

function displayPathComponents(records){
	var currK = ""
	var kHTML = "";
	var nHTML = "";
	var parent = document.getElementById('kegg-pathway-components');
	var componentDiv = document.createElement("div");
	componentDiv.className +='kegg-component';
	componentDiv.innerHTML = "<div class='component-kegg-data'><div class='ko-id'><b>KO Accession</b></div><div class='ko-description'><b>Description<b></div></div><div class='component-neurobase-data'><div class='component-neurobase-row'><div class='kegg-sb-cell'><b>NeuroBase ID</b></div><div class='kegg-annot-cell'><b>Sequence Annotation (E-value)</b></div></div>";
	parent.appendChild(componentDiv);
	//parse the records, displaying each
	for(i = 0; i < records.length; i++){
		r = records[i];
		//check if we're dealing with a new K0 term
		if(r['kegg_id'] != currK){
			if(currK != ""){
				nHTML += "</div>";
				//add the previous KO entry to the display list
				var componentDiv = document.createElement("div");
				componentDiv.className +='kegg-component';
				componentDiv.innerHTML = kHTML + nHTML;
				parent.appendChild(componentDiv);
			}
			//starting a new K0 entry in the display list
			kHTML = "<div class='component-kegg-data'><div class='ko-id'><a href='http://www.genome.jp/" + r['kegg_link'] + "' target='_blank'>" + r['kegg_long_id'] + "</a></div><div class='ko-description'>" + r['kegg_pathway_description'] + "</div></div>";
			nHTML = "<div class='component-neurobase-data'>";
			currK = r['kegg_id'];	
		}
		//add the neurobase data to the container
		nHTML += "<div class='component-neurobase-row'><div class='kegg-sb-cell'><a href='seqDetail.php?projectID=" + r['project_id'] + "&sbid=" + r['sb_id'] + "' target='_blank'>" + r['sb_id'] + "</a></div><div class='kegg-annot-cell'>" + r['annot'] + " (" + r['eval'] + ")</div></div>";	
	}
}