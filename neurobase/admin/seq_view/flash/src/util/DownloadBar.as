package util{
	
	import flash.display.Sprite;
	
	public class DownloadBar extends Sprite{
		
		private var bg:Sprite;
		private var bar:Sprite;
		private var maxWidth;
		
		public function DownloadBar(w:Number, h:Number, bgC:uint, barC:uint):void{
			this.maxWidth = w;
			bg = new Sprite();
			bar = new Sprite();
			bg.graphics.lineStyle(1, bgC, 1);
			bg.graphics.drawRect(0, 0, w, h);
			bar.graphics.beginFill(barC, 1);
			bar.graphics.drawRect(0, 0, w, h);
			bar.graphics.endFill();
			bar.width = 0;
			this.addChild(bar);
			this.addChild(bg);
		}
		
		public function setPercent(percent:Number):void{
			bar.width = Math.round(maxWidth * 100);
		}
		
	}
	
	
}