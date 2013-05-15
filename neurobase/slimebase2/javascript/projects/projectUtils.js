// JavaScript Document
var searchBox;
var selectedPfamPage;

window.onload = function(){
	
	var logo = document.getElementById('neurobase-logo');
	var leftMargin = (screen.width - 1000)/2; 
	logo.style.marginLeft = leftMargin + "px";	
	var navigation = document.getElementById('neurobase-nav');
	navigation.style.marginLeft = leftMargin + 624 + 12 + 30 + "px";
	var selectBox = document.getElementById('projectSelectElement');
	selectBox.onchange = loadProjectDetails;
	
	if (document.getElementsByClassName == undefined) {
		document.getElementsByClassName = function(className)
		{
			var hasClassName = new RegExp("(?:^|\\s)" + className + "(?:$|\\s)");
			var allElements = document.getElementsByTagName("*");
			var results = [];
			var element;
			for (var i = 0; (element = allElements[i]) != null; i++) {
				var elementClass = element.className;
				if (elementClass && elementClass.indexOf(className) != -1 && hasClassName.test(elementClass))
					results.push(element);
			}

			return results;
		}
	}
}

function loadProjectDetails(){
	var request = createAjaxRequest();
	var detailBox = $('projDetailContainer');
	var detailContent = document.getElementById('pSContentParagraph');
	var selectBox = document.getElementById('projectSelectElement');
	var index = selectBox.selectedIndex;
	var information = document.getElementById('main-information-content');
	information.style.display = 'none';
	currentSelected = "";
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var selectedID = selectBox.options[index].value;
		        var returnArray = JSON.parse(request.responseText);
			detailContent.innerHTML = returnArray[0];
			var organismContent = document.getElementById('organism-description');
			organismContent.innerHTML = returnArray[2];
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
			//clear the extra content container
			var contentContainer = document.getElementById('extraContentMain');
			contentContainer.innerHTML = "";
			//display the appropriate navigation options 
			displayNavigation(selectedID);
			//display the extra content
			displayExtraContent(selectedID);
			//change the background image
			var detailImg = document.getElementById('projectImage');
			detailImg.style.backgroundImage = 'url(../slimebase2/images/organisms/' + returnArray[1] + ')';
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadProjectDetails&args=" + selectBox.options[index].value, true);
	request.send(null);
}

function displayNavigation(selectedID){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var navContainer = document.getElementById('extraContentNav');
			navContainer.innerHTML = request.responseText;
			var navElements = document.getElementsByClassName('extraContentNavElement');
			for(i = 0; i < navElements.length; i++){
				if(navElements[i].id == 'go-link' || navElements[i].id == 'kegg-link'){
					continue;
				}
				navElements[i].onclick = function(){
					setActiveFeature(this.id);
					return(false);
				}
			}
			//display the project sequences
			initSequences(selectedID, 1, "none", "none");
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=displayNavigation&args=" + selectedID, true);
	request.send(null);
}	

//this function is called on page load
function initSequences(selectedID, pageNum, filter, filterValue){
	var request = createAjaxRequest();
	var contentContainer = document.getElementById('extraContentMain');
	contentContainer.innerHTML = "<div id='sequence-container'><div class='cross-comparison-header'>This tab displays sequences from a selected tissue/species sample.</div><div id='sequence-main'></div></div>";
	//select the sequence tab
	setActiveFeature("sequence-link");
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var sequenceContent = document.getElementById('sequence-main');
			sequenceContent.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/project_view/displaySequences.php?project_id=" + selectedID + "&pageNum=" + pageNum + "&filter=" + filter + "&filterValue=" + filterValue + "&sortVal=1", true);
	request.send(null);
}

//this function is called whenever sequence selections change
function displaySequences(selectedID, pageNum, filter, filterValue, sortVal){
	var request = createAjaxRequest();
	//select the sequence tab
	var navigationDisplays = document.getElementsByClassName('pageInfoContainer');
	for(i = 0; i < navigationDisplays.length; i++){
		navigationDisplays[i].innerHTML = "";	
	}
	var sequenceContent = document.getElementById('seqDisplay');
	sequenceContent.innerHTML = "";
	sequenceContent = document.getElementById('seq-header-container');
	sequenceContent.innerHTML = "";
	sequenceContent = document.getElementById('bottom-page-nav');
	sequenceContent.innerHTML = "";
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			sequenceContent = document.getElementById('sequence-main');
			sequenceContent.innerHTML = request.responseText;
		}
	}
	request.open("GET", "ajax/project_view/displaySequences.php?project_id=" + selectedID + "&pageNum=" + pageNum + "&filter=" + filter + "&filterValue=" + filterValue + "&sortVal=" + sortVal, true);
	request.send(null);
}

function displayExtraContent(projectID){
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			//parse the json encoding of the array and set the navigation 	
		    
			var returnArray = eval('(' + request.responseText + ')');
			var contentContainer = document.getElementById('extraContentMain');
			if(returnArray[0] != ""){
				contentContainer.style.minHeight = '700px';
				contentContainer.innerHTML += "<div id='cross-container'><div class='cross-comparison-header'>This cross-comparison tool allows you to view homologous or tissue-specific sequences from any two samples.</div><div id='cross-comparison-main'>" + returnArray[0] + "<div id='cross-comparison-result-container'></div></div></div>";
				
				//set the appropriately selected option for the first select box
				//var selectBox = document.getElementById('comparison-select-1');
				//switch the selection
				//selectBox.selectedIndex = returnArray[5];
				fetchCrossComparison(1);
			}
			if(returnArray[1] != ""){
				contentContainer.innerHTML += returnArray[1];
			}
			if(returnArray[2] != ""){
				contentContainer.innerHTML += returnArray[2];	
				searchBox = document.getElementById('pfam-search-element');
				searchBox.setAttribute("class", "idle-search");
				searchBox.onfocus = removeSearchGlass;
				searchBox.onkeyup = returnsCheckSearchText(projectID);
				searchBox.onblur = returnsCheckSearchText(projectID);
				searchBox.onkeydown = removeSearchGlass;
				selectedPfamPage = "A";
				attachPfamRowHandlers(projectID);
			}
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=loadExtraContent&args=" + projectID, true);
	request.send(null);
}


function checkSearchText(){
	if(this.value == ''){
		searchBox.setAttribute("class", "idle-search");	
	}
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			
		}
	}
	request.open("GET", "ajax/utils/pfam.php?method=updateSearchRes&args=" + projectID + "," + text, true);
	request.send(null);
}

function loadContent(id){
	var newActive = document.getElementById(id);
	var oldActive = document.getElementById(currentSelected);
	oldActive.style.backgroundImage = "";
	newActive.style.backgroundImage = "url(images/browse/extraSelected.jpg)";
	currentSelected = id;
}

function fetchCrossComparison(pageNum){
	var request = createAjaxRequest();
	var firstSelect = document.getElementById('comparison-select-1');
	var secondSelect = document.getElementById('comparison-select-2');
	if(firstSelect.length != 0){
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
	else{
		var resultContainer = document.getElementById('cross-comparison-result-container');
		resultContainer.innerHTML = "<br/><br/><h1>No cross-comparisons currently available, check back later.</h1>";	
	}
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