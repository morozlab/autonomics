package util{
	
	import flash.display.MovieClip;
	import flash.display.DisplayObject;
	import flash.events.MouseEvent;
	import flash.geom.Rectangle;
	
	public class VScrollBar extends MovieClip{
		
		var scrollBack:MovieClip = new MovieClip();
		var scrollPane:MovieClip = new MovieClip();
		var totalH:Number;
		var maskHeight:Number;
		var totalHeight:Number;
		var curScrollPosition:Number = 0;
		var scrollRatio:Number;
		var scrollT:DisplayObject;
		
		public function VScrollBar(target:DisplayObject, totalH:Number, maskHeight:Number){
			scrollBack.graphics.beginFill(0xE9E9E9, 1);
			scrollBack.graphics.drawRect(0, 0, 15, maskHeight);
			scrollBack.graphics.endFill();
			this.addChild(scrollBack);
			scrollPane.graphics.beginFill(0x9F9F9F, 1);
			scrollPane.graphics.drawRect(0, 0, 15, 40);
			scrollPane.graphics.endFill();
			scrollPane.name = 'scrollPane';
			this.addChild(scrollPane);
			scrollPane.addEventListener(MouseEvent.MOUSE_DOWN, startScroll);
			scrollPane.addEventListener(MouseEvent.MOUSE_UP, endScroll);
			this.totalH = totalH;
			this.maskHeight = maskHeight;
			this.initScroll();
			this.scrollT = target;
		}
		
		public function initScroll(){
			scrollPane.y = 0;
			curScrollPosition = 0;
			scrollRatio = (totalH - maskHeight)/(scrollBack.height - 40);
		}
		
		public function resizeScroll(newHeight:Number):void{
			this.totalH = totalH;
			this.initScroll();			
		}
		
		private function startScroll(e:MouseEvent){
			scrollPane.startDrag(false, new Rectangle(0, 0, 0, scrollBack.height - 40));
			this.stage.addEventListener(MouseEvent.MOUSE_UP, endScroll);
			this.stage.addEventListener(MouseEvent.MOUSE_MOVE, doScroll);
		}
		
		private function endScroll(e:MouseEvent){
			scrollPane.stopDrag();
			this.stage.removeEventListener(MouseEvent.MOUSE_UP, endScroll);
			this.stage.removeEventListener(MouseEvent.MOUSE_MOVE, doScroll);
		}

		/*internal function resizeScroll(opt:String, exp:Number){
			var tmp = MovieClip(scrollBar.getChildByName('scrollPane'));
			//determine the new scroll bar ratio
			scrollRatio = (totalH - maskHeight)/(maskHeight- 40);
			//figure out where the sequence box is in relation to the mask
			var dy:Number;
			if(opt == 'expand'){
				dy = Math.abs(sequenceBox.y - maskClip.y);
			   if(tmp.y != 0){
				   tmp.y = dy / scrollRatio;
			   }
			}
			else{
				if(sequenceBox.y + exp > 30){
					sequenceBox.y = 30;
				}
				else{
					sequenceBox.y += exp;
				}
				dy = Math.abs(sequenceBox.y - maskClip.y);
				if((dy / scrollRatio) > (scrollBar.height - tmp.height)){
					tmp.y = scrollBar.height - tmp.height;
				}
				else{
					tmp.y = dy / scrollRatio;
				}
			}
			//trace(maskClip.y + "|" + sequenceBox.y);
			curScrollPosition = tmp.y;
		}*/
		
		private function doScroll(e:MouseEvent){
			var posChange = curScrollPosition - scrollPane.y;
			curScrollPosition = scrollPane.y;
			scrollT.y = scrollT.y + (posChange * scrollRatio);
		}
		
		
	}
	
	
}