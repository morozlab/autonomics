package blast{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFieldType;
	import flash.events.MouseEvent;
	import sequencepanel.*;
	import blast.*;
	
	public class PastedSequence extends MovieClip{
		
		public var thisIndex:Number;
		
		public function PastedSequence(number:Number, ID:String, w:Number, bPanel:BlastPanel){
			
			var closeB:TextClip = new TextClip('x');
			closeB.x = 10;
			closeB.datum = this;
			closeB.addEventListener(MouseEvent.CLICK, (bPanel).removePastedSequence);
			var myText:TextField = new TextField();
			myText.x = closeB.x + closeB.width + 5;
			myText.text = ID;
			myText.autoSize = TextFieldAutoSize.LEFT;
			myText.width = w - closeB.width - closeB.x;
			myText.wordWrap = true;
			myText.multiline = true;
			this.addChild(myText);
			this.addChild(closeB);
			this.thisIndex = number;
			
		}
	}
	
	
}