package sequencepanel{
	
	import flash.display.Sprite;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
		
	public class DragSequence extends Sprite{
		
		public function DragSequence(myText:String){
			
			var tmpText:TextField = new TextField();
			tmpText.autoSize = TextFieldAutoSize.LEFT;
			tmpText.x = 0;
			tmpText.y = 0;
			tmpText.text = myText;
			this.addChild(tmpText);
			this.buttonMode = true;
			this.mouseChildren = false;
			
		}
		
	}
	
	
}