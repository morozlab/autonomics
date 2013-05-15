package blast{
	
	import flash.display.Sprite;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	
	public class QueryZone extends Sprite{
		
		private var queryList:Array = new Array();
		private var displayText:TextField = new TextField();
		private var displayBox:Sprite;
		private var head:BlastNode;
		private var tail:BlastNode;
		public var queryWidth;
		
		public function QueryZone(w:Number){
			
			this.name = 'query_zone';
			queryWidth = w;
			//create the border around zone
			displayBox = new Sprite();
			displayBox.graphics.lineStyle(1, 0x000000, 1);
			displayBox.graphics.drawRect(0, 0, w, w);
			this.addChild(displayBox);
			var tmp = new BlastNode("Drag a project here to set as query sequences", 0);
			this.addChild(tmp);
		}
		
		public function appendQuery(nodeName:String, nodeID:Number):void{
			if(this.getChildByName("Drag a project here to set as query sequences") != null){
				var tmp = this.getChildByName("Drag a project here to set as query sequences");
				this.removeChild(tmp);
			}
			var query = new BlastNode(nodeName, nodeID);
			if(this.head == null){
				this.head = query;
				this.tail = query;
				query.y = 0;
			}
			else{
				this.tail.nextNode = query;
				query.prevNode = this.tail;
				this.tail = query;
				query.y = query.prevNode.y + query.prevNode.height + 2;
			}
			this.addChild(query);
			
		}
		
		public function removeQuery(nodeID:Number):void{
			
		}
		
	}
	
}