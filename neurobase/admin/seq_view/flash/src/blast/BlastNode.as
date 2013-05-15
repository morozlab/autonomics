package blast{
	
	import flash.display.Sprite;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	
	public class BlastNode extends Sprite{
	
		private var blastText:TextField = new TextField();
		public var nextNode:BlastNode;
		public var prevNode:BlastNode;
		private var projID:Number;
		
		public function BlastNode(myText:String, myID:Number):void{
			
			this.mouseChildren = false;
			this.name = myText;
			
			blastText.text = myText;
			blastText.wordWrap = true;
			this.addChild(blastText);
			
		}
	}
}