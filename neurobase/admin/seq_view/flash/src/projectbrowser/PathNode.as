package projectbrowser{
	
	import flash.display.Sprite;
	import flash.display.MovieClip;
	import flash.text.*
	import flash.events.MouseEvent;
	import flash.net.Responder;
	import flash.net.NetConnection;
	import util.*;
	
	class PathNode extends Sprite{
		
		internal var nextNode;
		internal var prevNode;
		internal var textDisplay:TextField;
		internal var projID;
		internal var thisFormat:TextFormat = new TextFormat();

		public function PathNode(myText:String, myID:Number, clickable:Boolean):void{
			this.buttonMode = true;
			this.mouseChildren = false;
			this.name = myText;
			this.projID = myID;
			var tFont:TahomaFont = new TahomaFont();
			thisFormat.font = tFont.fontName;
			thisFormat.size = 14;
			textDisplay = new TextField();
			textDisplay.defaultTextFormat = thisFormat;
			textDisplay.embedFonts = false;
			textDisplay.htmlText = myText + "/";
			textDisplay.autoSize = TextFieldAutoSize.LEFT;
			this.addChild(textDisplay);
			if(clickable == true){
				enableClick();
			}
			this.addEventListener(MouseEvent.MOUSE_OVER, mouseOverText);
			this.addEventListener(MouseEvent.MOUSE_OUT, mouseOutText);
		}
		
		private function mouseOverText(e:MouseEvent):void{
			thisFormat.color = 0xF5B12B;
			textDisplay.setTextFormat(thisFormat);
			textDisplay.htmlText = "<u>" + textDisplay.htmlText + "</u>";
		}
		
		private function mouseOutText(e:MouseEvent):void{
			thisFormat.color = 0x000000;
			textDisplay.setTextFormat(thisFormat);			
			textDisplay.htmlText = this.name + "/";
		}
		
		protected function enableClick():void{
			this.addEventListener(MouseEvent.CLICK, doEnterProj);
		}
		
		protected function disableClick():void{
			this.removeEventListener(MouseEvent.CLICK, doEnterProj);
		}
		
		private function doEnterProj(e:MouseEvent):void{
			var path = PathDisplay(this.parent);
			var displayContainer = ProjectDisplay(path.parent);
			path.removeToNode(this.projID);
			displayContainer.enterProject(this.projID, this.name);
			//set this as the default project detail to display
			if(this.projID != -1){
				this.getProjDetails();
			}
		}
		
		private function getProjDetails():void{
			var responder:Responder = new Responder(setDefault);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query = "SELECT browser_img, browser_description FROM project_directory WHERE projectID ='" + this.projID + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function setDefault(res:Array):void{
			var displayContainer = ProjectDisplay(this.parent.parent);
			displayContainer.setFolderDetails(String(res[0]['browser_img']), String(res[0]['browser_description']), this.projID);
		}
	}
}