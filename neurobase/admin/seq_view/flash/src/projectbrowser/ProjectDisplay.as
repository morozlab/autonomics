package projectbrowser{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.events.*;
	// required to send/recieve data over AMF
	import flash.net.NetConnection;
	import flash.net.Responder;
	import projectbrowser.downloads.*;
	import util.*;
	import fl.transitions.Tween;
	import fl.transitions.easing.*;
	
	public class ProjectDisplay extends MovieClip{
		
		public var displayWidth = 300;
		private var connection:NetConnection;
		private var responder:Responder;
		private var query:String;
		private var titleText:TextField;
		private var numPermanent:Number = 0;
		private var searchBar:SearchBar;
		private var dlZone:DownloadZone;
		private var curFolderHeight = 0;
		private var oldFolder;
		internal var curFolder;
		public var quickBlastProject;
		public var quickBlastProjName;
		private var hasFiles:Boolean = false;
		private var curFolderName:String = "";
		public var path:PathDisplay;
		public var addToPath = "root";
		private var leftPadding:Number = 10;
		private var rightPadding:Number = 10;
		private var optionSpace:Number = 2;
		public var footer:DataSetBottom = new DataSetBottom();
		private var searchIcon:SearchIcon = new SearchIcon();
		private var dataOptions:MovieClip = new MovieClip();
		private var maskTween:Tween;
		private var cFasta:CustomFasta;
		private var projectDetail:ProjectDetail;
		
		//constructor of ProjectDisplay (duh)
		public function ProjectDisplay():void{
			//add the title text
			var dataHeader:DataSetHeader = new DataSetHeader();
			addChild(dataHeader);
			numPermanent++;
			//add the sequence options
			dataOptions.x = 15;
			dataOptions.y = dataHeader.y + dataHeader.height + 3;
			//add home button
			var homeButton:HomeButton = new HomeButton();
			homeButton.buttonMode = true;
			homeButton.addEventListener(MouseEvent.CLICK, goHome);
			dataOptions.addChild(homeButton);
			//dataOptions.addChild(downloadArrow);
			searchIcon.x = homeButton.x + homeButton.width + optionSpace - 3;
			searchIcon.addEventListener(MouseEvent.CLICK, displaySearchBox);
			searchIcon.buttonMode = true;
			searchIcon.y = -2;
			dataOptions.addChild(searchIcon);
			addChild(dataOptions);
			var tText:TextField = new TextField();
			numPermanent++;
			//instantiate the path display
			path = new PathDisplay(250);
			numPermanent++;
			path.y = dataOptions.y + dataOptions.height + 5;
			path.x = 17;
			addChild(path);
			projectDetail = new ProjectDetail();
			projectDetail.y = footer.height + 20;
			footer.addChild(projectDetail);
			footer.y = 130;
			addChild(footer);
			numPermanent++;
			//set the current display height
			resetFolderHeight();
			//addChild(footer);
			//begin displaying subprojects at top level
			
		}		
		
		/*private function setQuickBlast(e:MouseEvent):void{
			if(hasFiles == true){
				quickBlastProject = curFolder;
				quickBlastProjName = curFolderName;
				qbTween = new Tween(qbArea, 'y', Regular.easeOut, 0, (qbArea.height * -1), .5, true);
			}
			else{
				
			}
		}*/
		
		private function printY(e:Event):void{
			trace((MovieClip(footer.getChildByName('fastaCover'))).y);
		}
		
		private function displaySearchBox(e:MouseEvent):void{
			searchBar = new SearchBar();
			searchBar.x = dataOptions.x + dataOptions.width + 3;
			searchBar.y = dataOptions.y + 2;
			this.addChildAt(searchBar, 1);
			searchBar.startTween();
						
		}
		
		internal function enableSearchClick():void{
			searchIcon.removeEventListener(MouseEvent.CLICK, displaySearchBox);
			searchIcon.addEventListener(MouseEvent.CLICK, doSearchClick);
		}
		
		internal function disableSearchClick():void{
			searchIcon.addEventListener(MouseEvent.CLICK, displaySearchBox);
			searchIcon.removeEventListener(MouseEvent.CLICK, doSearchClick);
		}
		
		private function doSearchClick(e:MouseEvent):void{
			searchBar.hideBoxAndSearch();
		}
		
		private function goHome(e:MouseEvent){
			if(curFolder != 1){
				path.clearChildren();
				this.enterProject(1, 'projects');
		
			}
		}
		
		//handles the call to retrieve a list of projects from DB
		//translates a project name into a project ID
		/*public function enterProject(projName:String):void{
			responder = new Responder(handleProjEnter);
			connection = new NetConnection;
			connection.connect(inc.gateway);
			query = "SELECT projectID FROM project_directory WHERE project_name = '" + projName + "'";
			connection.call("DBI.execSelect", responder, query);
		}*/
		
		public function execSearch(search:String):void{
			//add that we're looking at search results
			if(path.getChildByName('Search Results') == null){
				var str = "Search Results";
				path.appendNode(str, -1, false);
				curFolder = -1;
			}
			responder = new Responder(displaySearchResult);
			connection = new NetConnection;
			connection.connect(MovieClip(root).gateway);
			query = "SELECT project_name, projectID, child_of, browser_img, browser_description FROM project_directory WHERE project_name LIKE '%" + search + "%' AND project_name != 'root'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function beginDL(e:MouseEvent):void{
			if(curFolder != -1){
				this.initDownload(this.curFolder, 'project');
			}
		}
		
		//get the path for the supplied project
		public function initDownload(projID:Number, type):void{
			responder = new Responder(compressFile);
			connection = new NetConnection;
			connection.connect(MovieClip(root).gateway);
			type = 'project';
			if(type == "FILE"){
				query = "SELECT path FROM project_files WHERE fileID = '" + projID + "'";
			}
			else{
				query = "SELECT path FROM project_directory WHERE projectID = '" + projID + "'";
			}
			connection.call("DBI.execSelect", responder, query);
		}
		
		//compresses file at the path returned from initDownload
		private function compressFile(res:Array):void{
			var path = String(res[0]['path']);
			var dl = DownloadZone(footer.getChildByName('dlZone'));
			dl.dlPath = path;
			responder = new Responder(dl.presentFile);
			connection = new NetConnection;
			connection.connect(MovieClip(root).gateway);
			connection.call("FlashDownload.compressFile", responder, path);
			dl.statText.text = "Compressing file for download, please wait.";
		}
		
		//uses retreived projectID to get a list of children, calls initDisplay to display the children
		public function enterProject(projID:Number, projName:String):void{
			if((projID != curFolder)){
				trace('calling enter project');
				responder = new Responder(changeDisplay);
				connection = new NetConnection;
				connection.connect(MovieClip(root).gateway);
				query = "SELECT project_name, path, child_of, projectID, browser_img, browser_description FROM project_directory WHERE  child_of ='" + projID + "'";
				connection.call("DBI.execSelect", responder, query);
				oldFolder = curFolder;
				curFolder = projID;
				curFolderName = projName;
			}
		}
		
		private function changeDisplay(res:Array):void{
			//clear the existing projects/files from the folder display
			clearChildren();
			//check that there was a result from DB
			if(res.length > 0 && oldFolder != -1){
				//clear the folder view
				//append new folder to existing path
				if(path.tailNode == null){
					path.appendNode(curFolderName, curFolder, true);				
				}
				else if(curFolder != path.tailNode.projID){
					path.appendNode(curFolderName, curFolder, true);
				}
				//add the children of the project
			}
			else{
				//this is a leaf in the project directory tree, so no projects to display
				//check if this we are entering from a search result
				if(oldFolder == -1){
					//clear the old path
					path.clearChildren();
					path.displayPath(curFolder);
				}
				else{
					//don't need to clear the path, just need to get the name of the newly entered project
					//and append it to the end
					path.appendNode(curFolderName, curFolder, true);
				}
			}
			var addedProject = 0;
			curFolderHeight = path.y + path.height + 5;
			for(var i = 0; i < res.length; i++){
					addedProject = 1;
					this.addNode(String(res[i]['project_name']), Number(res[i]['child_of']), Number(res[i]['projectID']), String(res[i]['browser_img']), String(res[i]['browser_description']));
			}
			if(addedProject == 0){
				//reset the folder height
				curFolderHeight = path.y + path.height + 5;
			}
			//display any files for the project
			displayProjectFiles(curFolder);
			//if this is the first project, add the download zone
			if(numPermanent == 3){
				dlZone = new DownloadZone();
				dlZone.name = 'dlZone';
				dlZone.x = 0;
				dlZone.y = path.stage.stageHeight - dlZone.height;
				this.addChild(dlZone);
				numPermanent++;
			}	
			
		}
		
		private function displayProjectFiles(projID:Number){
			responder = new Responder(doFileDisplay);
			connection = new NetConnection;
			connection.connect(MovieClip(root).gateway);
			query = "SELECT * FROM project_files WHERE projectID = '" + projID + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function doFileDisplay(res:Array):void{
			if(res.length > 0){
				hasFiles = true;
				var tmp = new LoadNode(String(res[0]['file_name']), String(res[0]['path']), Number(res[0]['projectID']));
				tmp.x = 20 + leftPadding;
				tmp.y = setNodeY(tmp);
				this.addChildAt(tmp, 1);
			}
				
		}
		
		private function displaySearchResult(res:Array):void{
			var tmp;
			clearChildren();
			//check if there was at least one search result from DB
			if(res.length > 0){
				for(var i = 0; i < res.length; i++){
					this.addNode(String(res[i]['project_name']), Number(res[i]['child_of']), Number(res[i]['projectID']), String(res[i]['browser_img']), String(res[i]['browser_description']));
				}
			}
			//if not, display that there were no matching results
			else{
				tmp = new ProjectNode('No Matching Results.', -1, -1, "none", "none");
				tmp.x = 20 + leftPadding;
				tmp.y = setNodeY(tmp);
				this.addChildAt(tmp, 1);
			}
			
		}
		
		private function addNode(projName:String, parentID:Number, projID:Number, img:String, desc:String):void{
			var tmp = new ProjectNode(projName, parentID, projID, img, desc);
			tmp.x = 20 + leftPadding;
			tmp.y = setNodeY(tmp);
			setFooter(tmp.y + tmp.height + 20);
			this.addChildAt(tmp, 1);
		}
		
		private function setFooter(amnt:Number):void{
			var footTween:Tween;
			if(amnt > 200){
				footer.y = amnt;
				//footTween = new Tween(footer, 'y', Regular.easeOut, footer.y, amnt, .5, true);
			}
			else{
				footer.y = 200;
				//footTween = new Tween(footer, 'y', Regular.easeOut, footer.y, 200, .5, true);
			}
		}
		
		private function clearChildren():void{
			while(this.numChildren > numPermanent){
				this.removeChildAt(1);
			}
			resetFolderHeight();
		}
		
		internal function resetFolderHeight(){
			curFolderHeight = path.y + path.height + 5;
		}
		
		private function setNodeY(node:ProjectNode):Number{
			var ret = curFolderHeight + 2;
			curFolderHeight = ret + node.height;
			return ret;
		}
		
		
		private function resizeDisplay():void{
			this.width = displayWidth;
		}
		
		//********************************************* Project Details Display Functions **************************//
		
		public function setFolderDetails(path:String, desc:String, projID:Number):void{
			projectDetail.setDefault(path, desc, projID);
		}
		
		public function displayProjectDetails(type:String, path:String, desc:String, projID:Number):void{
			if(type == "mouseOver"){
				projectDetail.loadAndShow(path, desc, projID);
			}
			else if(type == "mouseOut"){
				projectDetail.showDefault();
			}
		}		
		
	}
	
}

/*import flash.events.*;
// required to send/recieve data over AMF
import flash.net.NetConnection;
import flash.net.Responder;

var inc.gateway:String = "http://150.176.130.196:8888/amfphp/inc.gateway.php";
var connection:NetConnection;
var responder:Responder;
var query:String;
responder = new Responder(onResult, onFault);
connection = new NetConnection;
connection.connect(inc.gateway);

query = "SELECT * FROM project_directory";

connection.call("DBI.execSelect", responder, query);
function onResult(res:Array):void{
	result_text.text = res[0]['project_name'];
}

function onFault(fault:Object):void{
	result_text.text = 'error';
}*/