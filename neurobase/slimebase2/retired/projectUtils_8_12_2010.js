// JavaScript Document

function positionHeader(){
	var logo = document.getElementById('neurobase-logo');
	var leftMargin = (screen.width - 1000)/2; 
	logo.style.marginLeft = leftMargin + "px";	
	var navigation = document.getElementById('neurobase-nav');
	navigation.style.marginLeft = leftMargin + 624 + 12 + 30 + "px";
}


function loadProjectDetails(){
	var request = createAjaxRequest();
	var detailBox = $('projDetailContainer');
	var detailContent = document.getElementById('pSContentParagraph');
	var selectBox = document.getElementById('projectSelectElement');
	var index = selectBox.selectedIndex;
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var selectedID = selectBox.options[index].value;
			var returnArray = eval('(' + request.responseText + ')');
			detailContent.innerHTML = returnArray[0];
			var organismContent = document.getElementById('organism-description');
			organismContent.innerHTML = returnArray[2];
			//display the appropriate navigation options with their content
			displayExtraContent(selectedID);
			if(detailExpanded == false){
				detailExpanded = true;
				var boxTween = new Fx.Tween(detailBox);
				boxTween.set('display', 'block');
				boxTween.start('height', '0px', '250px');
				var extraContent = $('extraContent');
				var extraTween = new Fx.Tween(extraContent);
				extraTween.set('display', 'block');
				extraTween.start('opacity', '0', '1');
			}
			//change the background image
			var detailImg = document.getElementById('projectImage');
			detailImg.style.backgroundImage = 'url(../slimebase2/images/organisms/' + returnArray[1] + ')';
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadProjectDetails&args=" + selectBox.options[index].value, true);
	request.send(null);
}

function reloadProject(projectID){
	var request = createAjaxRequest();
	var detailBox = $('projDetailContainer');
	var detailContent = document.getElementById('pSContentParagraph');
	var selectBox = document.getElementById('projectSelectElement');
	var index = selectBox.selectedIndex;
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var selectedID = projectID;
			var returnArray = eval('(' + request.responseText + ')');
			if(detailExpanded == false){
				detailExpanded = true;
				var boxTween = new Fx.Tween(detailBox);
				boxTween.set('display', 'block');
				boxTween.start('height', '0px', '250px');
				var extraContent = $('extraContent');
				var extraTween = new Fx.Tween(extraContent);
				extraTween.set('display', 'block');
				extraTween.start('opacity', '0', '1');
			}
			detailContent.innerHTML = returnArray[0];
			//display the appropriate navigation options with their content
			displayExtraContent(selectedID);
			//change the background image
			var detailImg = document.getElementById('projectImage');
			detailImg.style.backgroundImage = 'url(../slimebase2/images/organisms/' + returnArray[1] + ')';
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadProjectDetails&args=" + projectID, true);
	request.send(null);	
}

window.onload = function(){
	document.title = "Cancer and Homarus PNAS Comparative Online Database";
}

function loadContent(id){
	var newActive = document.getElementById(id);
	var oldActive = document.getElementById(currentSelected);
	oldActive.style.backgroundImage = "";
	newActive.style.backgroundImage = "url(images/browse/extraSelected.jpg)";
	currentSelected = id;
}

function displayExtraContent(projectID){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			//parse the json encoding of the array and set the navigation 	
			var returnArray = eval('(' + request.responseText + ')');
			var navContainer = document.getElementById('extraContentNav');
			var contentContainer = document.getElementById('extraContentMain');
			navContainer.innerHTML = "";
			var toSelect = "";
			navContainer.innerHTML += returnArray[4];
			if(returnArray[0] != ""){
				navContainer.innerHTML += returnArray[0];
				contentContainer.style.minHeight = '500px';
				contentContainer.innerHTML = "<div id='cross-container'><div class='cross-comparison-header'>This cross-comparison tool allows you to view homologous or tissue-specific sequences from any two samples.</div><div id='cross-comparison-main'>" + returnArray[1] + "<div id='cross-comparison-result-container'></div></div></div>";
				//set the appropriately selected option for the first select box
				//var selectBox = document.getElementById('comparison-select-1');
				//switch the selection
				//selectBox.selectedIndex = returnArray[5];
				fetchCrossComparison(1);
				toSelect = 'cross-link';
			}
			else{
				contentContainer.innerHTML = "";	
			}
			if(returnArray[2] != ""){
				navContainer.innerHTML += returnArray[2];
				contentContainer.innerHTML += returnArray[3];
				if(toSelect == ""){
					toSelect = 'abstract-link';				
				}
			}
			if(toSelect != ""){
				setActiveFeature(toSelect);	
			}
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadExtraContent&args=" + projectID, true);
	request.send(null);
}

function fetchCrossComparison(pageNum){
	var request = createAjaxRequest();
	var firstSelect = document.getElementById('comparison-select-1');
	var secondSelect = document.getElementById('comparison-select-2');
	var selected1 = firstSelect.options[firstSelect.selectedIndex].value;
	var selected2 = secondSelect.options[secondSelect.selectedIndex].value;
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var resultContainer = document.getElementById('cross-comparison-result-container');
			resultContainer.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/project_view/populateComparisonTable.php?page=" + pageNum + "&project1=" + selected1 + "&project2=" +selected2, true);
	request.send(null);
}

function loadSelectedPage(){
	var pageSelect = document.getElementById('pageSelect');
	fetchCrossComparison(pageSelect.selectedIndex + 1);	
}

function loadComparisonAlignment(projectID1, projectID2, sbID1, sbID2){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var container = document.getElementById(projectID1 + "-" + projectID2 + "-" + sbID1 + "-" + sbID2);
			container.style.background = '#f3f3f3';
			container.style.border = '1px dashed #bcdae5';
			container.style.display = 'block';
			container.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadAlign&args=" + projectID1 + "," + projectID2 + "," + sbID1 + "," + sbID2, true);
	request.send(null);
}

function setActiveFeature(id){
	if(currentSelected != id && currentSelected != ""){
		var idArray = currentSelected.split("-");
		var linkToHide = document.getElementById(currentSelected);
		var contentToHide = document.getElementById(idArray[0] + '-container');
		linkToHide.style.backgroundImage = "";
		contentToHide.style.display = "none";
	}
	var activeLink = document.getElementById(id);
	var splitId = id.split("-");
	var activeContent = document.getElementById(splitId[0] + "-container");
	activeContent.style.display = 'block';
	activeLink.style.backgroundImage = "url(images/browse/extraSelected.jpg)";
	currentSelected = id;
}

function reloadSecondSelect(){
	var secondSelect = document.getElementById('comparison-select-2-span');
	var firstSelect = document.getElementById('comparison-select-1');
	var selectedVal = firstSelect.options[firstSelect.selectedIndex].value;
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			secondSelect.innerHTML = "<select name='comparison-select-2' id='comparison-select-2' class='comparison-select' onchange='fetchCrossComparison(1)'>" + request.responseText + "</select>";
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=reloadSecondSelect&args=" + selectedVal, true);
	request.send(null); 
}