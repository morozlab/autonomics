package sequencepanel{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	
	public class TextClip extends MovieClip{
		
		public var thisText:TextField = new TextField();
		public var datum;
		public function TextClip(myText:String){
			
			var thisFormat:TextFormat = new TextFormat();
			var thisFont:Font = new Font1();
			thisFormat.font = thisFont.fontName;
			thisFormat.size = 11;
			thisText.defaultTextFormat = thisFormat;
			thisText.embedFonts = true;
			thisText.antiAliasType = AntiAliasType.ADVANCED;
			thisText.x = 0;
			thisText.y = 0;
			thisText.text = myText;
			thisText.autoSize = TextFieldAutoSize.LEFT;
			this.mouseChildren = false;
			this.buttonMode = true;
			this.addChild(thisText);

		}
	}
	
}