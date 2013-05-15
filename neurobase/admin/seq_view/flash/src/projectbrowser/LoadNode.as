package projectbrowser{
		
	import flash.events.MouseEvent;
	import flash.events.Event;
	import flash.text.TextFieldAutoSize;
	import flash.display.Sprite;
	import flash.display.Loader;
	import flash.net.URLRequest;
	import sequencepanel.*;
	
	public class LoadNode extends ProjectNode{
		
		internal var loadPath:String;
		internal var fileIcon:Sprite = new Sprite();
		internal var img:Loader;
		internal var ID;
		
		public function LoadNode(myName:String, myPath:String, myID:Number){
			super(myName, -1, myID, "none", "none");
			img = new Loader();
			var req = new URLRequest("../images/file_icon.jpg");
		    img.contentLoaderInfo.addEventListener(Event.COMPLETE, attachIcon);
			img.load(req);
			this.loadPath = myPath;
			this.ID = myID;
			//left-justify the text
			projectText.autoSize = TextFieldAutoSize.LEFT;
			projectText.text = 'Load Sequences';
			this.addChild(projectText);
			this.buttonMode = true;
			this.mouseChildren = false;
			//this.doubleClickEnabled = true;
			this.addEventListener(MouseEvent.CLICK, this.doSeqLoad);
			this.removeEventListener(MouseEvent.ROLL_OVER, displayThisDetail);
			this.removeEventListener(MouseEvent.ROLL_OUT, displayFolderDetail);	
			this.removeEventListener(MouseEvent.CLICK, changeDisplay);
		}
		
		private function doSeqLoad(e:MouseEvent):void{
			var tmp:SequencePanel = SequencePanel(this.parent.parent.getChildByName('seqPanel'));
			tmp.loadSeqs(this.ID, 'project', 1);
		}
		
		protected function attachIcon(e:Event):void{
			this.projectText.x = img.x + img.width;
			this.addChild(img);
		}
		
		protected override function monitorMouse(e:MouseEvent):void{
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_UP, changeDisplay);
			e.currentTarget.addEventListener(MouseEvent.MOUSE_MOVE, doDrag);
			this.mouseTimer.start();
			e.currentTarget.addEventListener(Event.ENTER_FRAME, checkTimer);
		}
		
		protected override function selfDestruct(e:MouseEvent):void{
			//check if user dragged the project over the download area
			/*if(checkDL()){
				//user wants to download a project, initiate download
				var dlZone = parent.getChildByName('dlZone');
				if(dlZone.checkEnabled()){
					trace('moo');
					dlZone.disable();
					var rent = ProjectDisplay(this.parent);
					rent.initDownload(this.projID, "FILE");
				}
			}
			else if(nodeHitTest(this.parent.parent.getChildByName('blast_panel'))){
				
			}
			parent.stage.removeEventListener(Event.MOUSE_LEAVE, outOfBounds);
			parent.stage.removeEventListener(MouseEvent.MOUSE_UP, selfDestruct);
			parent.removeChild(FileNode(parent.getChildByName(this.name)));
			*/
		}
		
		
		public override function doClone():ProjectNode{
			var tmpNode = new LoadNode((this.name + 'clone'), this.loadPath, this.ID);
			var tmpImg = new Loader();
			var tmpURL = new URLRequest('images/file_icon.jpg');
			tmpImg.load(tmpURL);
			tmpNode.addChild(tmpImg);
			tmpNode.projectText.text = this.name;
			tmpNode.x = this.x;
			tmpNode.y = this.y;
			return tmpNode;
		}

	}
	
}