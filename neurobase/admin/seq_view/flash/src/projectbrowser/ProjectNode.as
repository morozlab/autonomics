package projectbrowser{
	
	import flash.display.MovieClip;
	import flash.display.DisplayObject;
	import flash.text.TextField;
	import flash.text.TextFormat;
	import flash.text.TextFieldAutoSize;
	import flash.text.AntiAliasType;
	import flash.events.MouseEvent;
	import flash.utils.Timer;
	import flash.events.Event;
	import sequencepanel.*;
	
	public class ProjectNode extends MovieClip{
		//keeps track of where this 
		var childOf:Number;
		//displays the actual text of the node
		var projectText:TextField = new TextField();
		var projID:Number;
		var mouseTimer:Timer = new Timer(100);
		var thisFormat:TextFormat = new TextFormat();
		public var desc:String;
		public var path:String;
		
		//initializes each project identifier
		public function ProjectNode(myName:String, myParent:Number, myID:Number, myPath:String, myDesc:String){
			this.name = myName;
			this.childOf = myParent;
			this.projID = myID;
			this.path = myPath;
			this.desc = myDesc;
			var projectFolder:ProjectFolder = new ProjectFolder();
			this.addChild(projectFolder);
			//left-justify the text
			var tFont:TahomaFont = new TahomaFont();
			thisFormat.font = tFont.fontName;
			thisFormat.size = 13;
			projectText.defaultTextFormat = thisFormat;
			projectText.embedFonts = false;
			projectText.autoSize = TextFieldAutoSize.LEFT;
			projectText.text = myName;
			projectText.antiAliasType = AntiAliasType.NORMAL;
			projectText.x = projectFolder.x + projectFolder.width + 3;
			this.addChild(projectText);
			this.buttonMode = true;
			this.mouseChildren = false;
			//this.doubleClickEnabled = true;
			this.addEventListener(MouseEvent.MOUSE_OVER, mouseOverText);
			this.addEventListener(MouseEvent.MOUSE_OUT, mouseOutText);
			//this.addEventListener(MouseEvent.MOUSE_DOWN, monitorMouse);
			//add event listeners for displaying the images for this node
			this.addEventListener(MouseEvent.ROLL_OVER, displayThisDetail);
			this.addEventListener(MouseEvent.ROLL_OUT, displayFolderDetail);
			this.addEventListener(MouseEvent.CLICK, changeDisplay);
		}
		
		protected function displayThisDetail(e:MouseEvent):void{
			ProjectDisplay(this.parent).displayProjectDetails('mouseOver', this.path, this.desc, this.projID);
		}
		
		protected function displayFolderDetail(e:MouseEvent):void{
			if(this.parent != null){
				ProjectDisplay(this.parent).displayProjectDetails('mouseOut', this.path, this.desc, this.projID);
			}
		}
		
		private function mouseOverText(e:MouseEvent):void{
			thisFormat.color = 0x0069CA;
			projectText.setTextFormat(thisFormat);
			
		}
		
		private function mouseOutText(e:MouseEvent):void{
			thisFormat.color = 0x000000;
			projectText.setTextFormat(thisFormat);
		}
		
		public function doClone():ProjectNode{
			var tmpNode = new ProjectNode((this.name + 'clone'), this.childOf, this.projID, this.path, this.desc);
			tmpNode.projectText.text = this.name;
			tmpNode.x = this.x;
			tmpNode.y = this.y;
			return tmpNode;
		}
		
		protected function doDrag(e:MouseEvent):void{
			e.currentTarget.removeEventListener(Event.ENTER_FRAME, checkTimer);
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_MOVE, doDrag);
			this.mouseTimer.reset();
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_UP, changeDisplay);
			this.cloneAndDrag();
		}
		
		protected function monitorMouse(e:MouseEvent):void{
			e.currentTarget.addEventListener(MouseEvent.MOUSE_UP, changeDisplay);
			e.currentTarget.addEventListener(MouseEvent.MOUSE_MOVE, doDrag);
			this.mouseTimer.start();
			e.currentTarget.addEventListener(Event.ENTER_FRAME, checkTimer);
		}
		
		protected function checkTimer(e:Event):void{
			if(this.mouseTimer.currentCount >= 4){
				e.currentTarget.removeEventListener(MouseEvent.MOUSE_UP, changeDisplay);
				e.currentTarget.removeEventListener(MouseEvent.MOUSE_MOVE, doDrag);
				this.mouseTimer.reset();
				this.cloneAndDrag();
			}
		}
		
		protected function changeDisplay(e:MouseEvent):void{
			//e.currentTarget.removeEventListener(Event.ENTER_FRAME, checkTimer);
			//e.currentTarget.removeEventListener(MouseEvent.MOUSE_MOVE, doDrag);
			ProjectDisplay(this.parent).setFolderDetails(this.path, this.desc, this.projID);
			this.mouseTimer.reset();
			var myParent = ProjectDisplay(this.parent);
			myParent.enterProject(this.projID, this.projectText.text);
		}
		
		protected function cloneAndDrag():void{
			var tmp:ProjectNode = this.doClone();
			parent.addChild(tmp);
			tmp.startDrag();
			tmp.addEventListener(MouseEvent.MOUSE_UP, tmp.selfDestruct);
			tmp.stage.addEventListener(MouseEvent.MOUSE_UP, tmp.selfDestruct);
			tmp.stage.addEventListener(Event.MOUSE_LEAVE, tmp.outOfBounds);
		}
		
		protected function outOfBounds(e:Event):void{
			parent.removeChild(parent.getChildByName(this.name));
			this.stage.removeEventListener(Event.MOUSE_LEAVE, outOfBounds);
		}
		
		protected function selfDestruct(e:MouseEvent):void{
			//check if user dragged the project over the download area
			if(checkDL()){
				//user wants to download a project, initiate download
				var dlZone = parent.getChildByName('dlZone');
				if(dlZone.checkEnabled()){
					dlZone.disable();
					var rent = ProjectDisplay(this.parent);
					rent.initDownload(this.projID, "PROJECT");
				}				
			}
			else if(parent.getChildByName(this.name).hitTestObject(parent.parent.getChildByName('seqPanel'))){
				//load the sequences for this project
				var seqPanel = SequencePanel(this.parent.parent.getChildByName('seqPanel'));
				seqPanel.loadSeqs(this.projID, "project", 1);
			}
			parent.stage.removeEventListener(Event.MOUSE_LEAVE, outOfBounds);
			parent.stage.removeEventListener(MouseEvent.MOUSE_UP, selfDestruct);
			parent.removeChild(ProjectNode(parent.getChildByName(this.name)));
		}
		
		
		protected function checkDL():Boolean{
			if(this.parent.getChildByName('dlZone')){
				var testObj = parent.getChildByName('dlZone');
				var cloneObj = parent.getChildByName(this.name);
				if(cloneObj.hitTestObject(testObj)){
					return true;
				}
				else{
					return false;
				}
			}
			else{
				return false;
			}
		}
		
		
	}
}