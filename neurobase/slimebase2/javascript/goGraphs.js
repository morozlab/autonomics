// JavaScript Document

var oldDisplay;

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

function gatherGraphData(){

	var request = createAjaxRequest();
	
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var topLevelArray = eval(request.responseText);
			var process = parseInt(topLevelArray[0]);
			var molec = parseInt(topLevelArray[1]);
			var cell = parseInt(topLevelArray[2]);
			var options = PlotKit.Base.officeCyan();
			options['xTicks'] = new Array();
			options['xTicks'] = [{v:0, label:"Biological Process"}, {v:1, label:"Mollecular Function"}, {v:2, label:"Cellular Component"}];
			options['padding'] = {left: 10, right: 10, top: 10, bottom: 30};
			options['axisLabelWidth'] = 70;
			options['axisLabelFontSize'] = 11;
			options['drawBackground'] = false;
			var layout = new PlotKit.Layout("pie", options);
			layout.addDataset("Top Level", [[0, process], [1, molec], [2, cell]]);
			layout.evaluate();
			var canvas = document.getElementById("top_categories_canvas");
    		var plotter = new PlotKit.SweetCanvasRenderer(canvas, layout, options);
    		plotter.render();
			//call function to graph other levels of the tree
			graphSlim("P", options);
			
		}
	}
	request.open("GET", '../ajax/gatherGraphData.php?projectID=' + projectID + '&level=top', true);
	request.send(null);

}

function graphSlim(category){
	
	var request = createAjaxRequest();
	
	var target;
	if(category == "P"){
		target = "bio_p";
	}
	if(category == "F"){
		target = "mollec";	
	}
	if(category == "C"){
		target = "cell";
	}
	
	//hide divs and headers
	var toHide = document.getElementById(target + "_div");
	var h1 = document.getElementById(target + "_h1");
	toHide.style.display = 'none';
	h1.style.display = 'none';
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var resultArray = eval('(' + request.responseText + ')');
			var tickArray = new Array();
			var dataArray = new Array();
			var index = 0;
			for(var key in resultArray){
				tickArray[index] = {v:index, label:key.substring(0, 40)};
				dataArray[index] = [index, parseInt(resultArray[key])];
				//alert([index, resultArray[key]]);
				index++;
				
			}
			//calculate % each bar is
			var options = PlotKit.Base.officeCyan();
			options['padding'] = {left: 50, right: 10, top: 10, bottom: 30};
			options['axisLabelWidth'] = 60;
			options['axisLabelFontSize'] = 11;
			options['drawBackground'] = false;
			//set the tick labels to the tickArray
			options['xTicks'] = tickArray;
			var layout = new PlotKit.Layout("bar", options);
			layout.addDataset(target, dataArray);
			layout.evaluate();
			
			var canvas = document.getElementById(target + "_canvas");
    		var plotter = new PlotKit.SweetCanvasRenderer(canvas, layout, options);
    		plotter.render();
			//call function to graph other levels of the tree
			if(category == "P"){				
				graphSlim("F");
			}
			else if(category == "F"){
				graphSlim("C");	
			}
		}
	}
	request.open("GET", '../ajax/gatherGraphData.php?projectID=' + projectID + '&level=' + category, true);
	request.send(null);
	
}

function toggleDisplay(name){
	var element = document.getElementById(name + "_div");
	var header = document.getElementById(name + "_h1");
	if(element.style.display != 'none'){
		element.style.display = 'none';
		header.style.display = 'none';
	}
	else{
		element.style.display = 'block';
		header.style.display = 'block';
	}
}

//functions to count and init the first highest-level GO categories
			function initTopLevel(){
				//create the ajax request
				var request = createAjaxRequest();				

				//function that will handle the response from the script that creates the top-level GO categories
				request.onreadystatechange = function(){
					if(request.readyState == 4){
						var topLevelDisplay = document.getElementById('top_level');
						topLevelDisplay.innerHTML = request.responseText;
					}
				}
				
				//send the ajax request
				request.open("GET", "../ajax/initTopLevel.php?projectID=" + projectID, true);
				request.send(null); 
			}
			
//function to expand or collapse the individual levels of the GO Tree
			function toggleView(target, mode, level, expected){
				var targetEl = document.getElementById(target);
				var targetDiv = document.getElementById(target + "_div");
				var childDiv = document.getElementById(target + "_children");
				if(mode == 'expand'){
					//let the user know we're loading the GO info
					var url = "../ajax/expandGOCategory.php?level=" + level + "&projectID=" + projectID + "&term=" + target + "&expected=" + expected;	
					//send the appropriate request to expand the top-level GO category
					var req = createAjaxRequest();	
					//handle the data retunred from the server
					req.onreadystatechange = function(){
						if(req.readyState == 4){
							//change the link to collapse mode
							targetDiv.innerHTML = "<a target='_self' href='#' id='" + target + "' onclick='toggleView(this.id, \"collapse\", \"" + level + "\");return(false);'>" + targetEl.innerHTML + "</a><div id='" + target + "_children'></div>";
							childDiv = document.getElementById(target + "_children");					
							childDiv.innerHTML = req.responseText;
						}	
					}
						
					//send request to expand the top-level
					req.open("GET", url , true);
					req.send(null);
					
					
				}
				else if(mode == 'collapse'){
					//remove the content in the child div and change the link
					childDiv.innerHTML = "";
					targetDiv.innerHTML = "<a href='#' id=" + target + " onclick='toggleView(this.id, \"expand\", \"" + level + "\");return(false);'>" + targetEl.innerHTML + "<a><div id='" + target + "_children></div>";
				}
			}

//function to expand or collapse the individual levels of the GO Tree
function toggleCategory(target, mode, extra){
	var targetEl = document.getElementById(target);
	var targetDiv = document.getElementById(target + "_div");
	var childDiv = document.getElementById(target + "_children");
	if(mode == 'expand'){
		//let the user know we're loading the GO info
		var url = "../ajax/expandGOCategory.php?extra=" + extra + "&projectID=" + projectID + "&term=" + target;	
		//send the appropriate request to expand the top-level GO category
		var req = createAjaxRequest();	
		//handle the data retunred from the server
		req.onreadystatechange = function(){
			if(req.readyState == 4){
			//change the link to collapse mode
				targetDiv.innerHTML = "<a target='_self' href='#' id='" + target + "' onclick='toggleCategory(this.id, \"collapse\", \"" + extra + "\");return(false);'>" + targetEl.innerHTML + "</a><div id='" + target + "_children'></div>";
				childDiv = document.getElementById(target + "_children");					
				childDiv.innerHTML = req.responseText;
			}	
		}
						
		//send request to expand the top-level
		req.open("GET", url , true);
		req.send(null);
	}
	else if(mode == 'collapse'){
	//remove the content in the child div and change the link
		childDiv.innerHTML = "";
		targetDiv.innerHTML = "<a href='#' id=" + target + " onclick='toggleCategory(this.id, \"expand\", \"" + extra + "\");return(false);'>" + targetEl.innerHTML + "<a><div id='" + target + "_children></div>";
	}
}

function resetChecks(){
		
}