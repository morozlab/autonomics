package projectbrowser{
	
	import flash.display.Sprite;
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.net.NetConnection;
	import flash.net.Responder;
	
	class PathDisplay extends Sprite{
		
		internal var headNode;
		internal var tailNode;
		private var maxWidth;
		private var curLineWidth = 0;
		var startingTitle:PathFolder;
		private var displayArray:Array;
		private var displayLength;
		private var responder:Responder;
		private var connection:NetConnection;
		private var baseWidth:Number;
		
		public function PathDisplay(swidth:Number):void{
			maxWidth = swidth;
			startingTitle = new PathFolder();
			startingTitle.visible = false;
			this.addChild(startingTitle);
			curLineWidth = 0;
			baseWidth = curLineWidth;
		}
		
		public function displayPath(parentID:Number):void{
			//build the path display array from the bottom up
			displayArray = new Array();
			responder = new Responder(buildProjectPath);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query = "SELECT project_name, projectID, child_of FROM project_directory WHERE projectID = '" + parentID + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function buildProjectPath(res:Array):void{
			var tmpArray:Array = new Array();
			//store the project name and ID in temporary array to add to display array
			tmpArray[0] = String(res[0]['project_name']);
			tmpArray[1] = Number(res[0]['projectID']);
			displayArray.unshift(tmpArray);
			//two cases, the first is that this is the projects directory, and we can stop
			//the other case is that it is some other directory, and we need to keep moving up
			if(String(res[0]['project_name']) == 'projects'){
				this.doDisplay();
			}
			else{
				responder = new Responder(buildProjectPath);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				var query = "SELECT project_name, projectID, child_of FROM project_directory WHERE projectID = '" + Number(res[0]['child_of']) + "'";
				connection.call("DBI.execSelect", responder, query);
			}
		}
		
		private function doDisplay():void{
			for(var i = 0; i < displayArray.length; i++){
				this.appendNode(displayArray[i][0], displayArray[i][1], true);
			}
			//set the proper location of the folder items
			this.parent.positionProjectOptions();
		}
		
		public function appendNode(pathName:String, ID:Number, clickable:Boolean):void{
			var tmp = new PathNode(pathName, ID, clickable);
			if(curLineWidth + tmp.width < maxWidth){
				if(tailNode != null){
					tmp.x = tailNode.x + tailNode.width;
					tmp.y = tailNode.y;
					curLineWidth += tmp.width;
					
				}
				else{
					tmp.x = 0;
					tmp.y = 0;
					curLineWidth += tmp.width;
				}
			}
			else{
					curLineWidth = tmp.width;
					tmp.x = 0; 
					tmp.y = this.height + 2;
			}
			trace(curLineWidth);
			tmp.x -= 1;
			ProjectDisplay(parent).resetFolderHeight();
			doAppend(tmp);
			this.addChildAt(tmp, 1);
		}
		
		public function clearChildren():void{
			while(this.numChildren > 1){
				this.removeChildAt(1);
			}
			this.curLineWidth = 0;
			headNode = null;
			tailNode = null;
		}
		
		public function removeToNode(nodeID:Number){
			while(tailNode.projID != nodeID){
				var tmpNode = tailNode;
				tailNode = tailNode.prevNode;
				this.removeChild(tmpNode);
				if(curLineWidth - tmpNode.width <= 0){
					curLineWidth = tailNode.x + tailNode.width;
				}
				else{
					curLineWidth -= tmpNode.width;
				}
				ProjectDisplay(parent).curFolder = -2;
			}
			
		}
				
		
		private function doAppend(tmp:PathNode):void{
			if(headNode == null || tailNode == null){
				headNode = tmp;
				tailNode = tmp;
			}
			else{
				tmp.prevNode = tailNode;
				tailNode.nextNode = tmp;
				tailNode = tmp;
			}
		}
	}
	
	
}