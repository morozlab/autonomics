package util{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.AntiAliasType;
	import flash.text.TextFormat;
	import flash.text.Font;
	import fl.transitions.Tween;
	import fl.transitions.TweenEvent;
	import fl.transitions.easing.*;
	
	public class LoadingMovie extends MovieClip{
		
		private var sideSize:Number;
		private var sideSpeed:Number = .75
		private var colorArray:Array = new Array(uint(0x36E600), uint(0xE600C8), uint(0x00E6DE), uint(0xE60020), uint(0xE69800), uint(0x8F8F8F), uint(0x000000));
		//create the four sides
		var left:MovieClip = new MovieClip();
		var right:MovieClip = new MovieClip();
		var top:MovieClip = new MovieClip();
		var bottom:MovieClip = new MovieClip();
		//array to hold all the tweens so they don't get out of hand
		var tweenArray:Array = new Array();
		//text field to hold the text thingy
		var workingText:TextField = new TextField();
		
		public function LoadingMovie(size:Number, txt:String):void{
			sideSize = size;
			var workFormat:TextFormat = new TextFormat();
			var workFont:Font = new Font1();
			workFormat.font = workFont.fontName;
			workFormat.size = 12;
			workFormat.color = 0x000000;
			//workFormat.bold = true;
			workingText.defaultTextFormat = workFormat;
			workingText.embedFonts = true;
			workingText.antiAliasType = AntiAliasType.ADVANCED;
			workingText.autoSize = TextFieldAutoSize.LEFT;
			workingText.x = sideSize + 5;
			workingText.text = txt;
			workingText.y = (sideSize - workingText.height)/2;
			this.addChild(workingText);
			//make everything invisible, cept for yon text field
			left.visible = false;
			right.visible = false;
			top.visible = false;
			bottom.visible = false;
			//add all the childrens
			this.addChild(left);
			this.addChild(right);
			this.addChild(top);
			this.addChild(bottom);
		}
		
		public function startAnimation():void{
			//clear the movie clips
			clearClips();
			//set the coordinates that won't change during tweening
			left.x = 0;
			right.x = sideSize;
			right.y = 0;
			top.y = 0;
			top.x = 0;
			
			bottom.y = sideSize;
			//pick a random color
			var randIndex:Number = Math.round(Math.random() * 6);
			var color:uint = colorArray[randIndex];
			//draw the  graphics
			left.visible = false;
			left.graphics.beginFill(color, 1);
			left.graphics.drawRect(0, 0, 1, this.sideSize);
			left.graphics.endFill();
			right.visible = false;
			right.graphics.beginFill(color, 1);
			right.graphics.drawRect(0, 0, 1, this.sideSize);
			right.graphics.endFill();
			top.visible = false;
			top.graphics.beginFill(color, 1);
			top.graphics.drawRect(0, 0, this.sideSize, 1);
			top.graphics.endFill();
			bottom.visible = false;
			bottom.graphics.beginFill(color, 1);
			bottom.graphics.drawRect(0, 0, this.sideSize, 1);
			bottom.graphics.endFill();
			//start the tween chain, whee			
			left.visible = true;
			tweenArray[0] = new Tween(left, 'y', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[1] = new Tween(left, 'height', Regular.easeOut, 1, sideSize, sideSpeed, true);
			tweenArray[1].addEventListener(TweenEvent.MOTION_FINISH, tweenTop);
		}
		
		public function setTextFormat(format:TextFormat):void{
			this.workingText.setTextFormat(format);
		}
		
		private function tweenTop(e:TweenEvent):void{
			tweenArray[1].removeEventListener(TweenEvent.MOTION_FINISH, tweenTop);
			top.visible = true;
			tweenArray[2] = new Tween(top, 'width', Regular.easeOut, 0, sideSize, sideSpeed, true);
			tweenArray[2].addEventListener(TweenEvent.MOTION_FINISH, tweenRight);
		}
		
		private function tweenRight(e:TweenEvent):void{
			tweenArray[2].removeEventListener(TweenEvent.MOTION_FINISH, tweenRight);
			right.visible = true;
			tweenArray[3] = new Tween(right, 'height', Regular.easeOut, 0, sideSize, sideSpeed, true);
			tweenArray[3].addEventListener(TweenEvent.MOTION_FINISH, tweenBottom);
		}
		
		private function tweenBottom(e:TweenEvent):void{
			tweenArray[3].removeEventListener(TweenEvent.MOTION_FINISH, tweenBottom);
			bottom.visible = true;
			tweenArray[4] = new Tween(bottom, 'x', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[5] = new Tween(bottom, 'width', Regular.easeOut, 0, sideSize + 1, sideSpeed, true);
			tweenArray[5].addEventListener(TweenEvent.MOTION_FINISH, removeLeft);
		}
		
		private function removeLeft(e:TweenEvent):void{
			tweenArray[5].removeEventListener(TweenEvent.MOTION_FINISH, removeLeft);
			tweenArray[1] = new Tween(left, 'height', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[1].addEventListener(TweenEvent.MOTION_FINISH, removeTop);
		}
		
		private function removeTop(e:TweenEvent):void{
			tweenArray[1].removeEventListener(TweenEvent.MOTION_FINISH, removeTop);
			left.visible = false;
			tweenArray[2] = new Tween(top, 'x', Regular.easeOut, 0, sideSize, sideSpeed, true);
			tweenArray[3] = new Tween(top, 'width', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[3].addEventListener(TweenEvent.MOTION_FINISH, removeRight);
		}
		
		private function removeRight(e:TweenEvent):void{
			tweenArray[3].removeEventListener(TweenEvent.MOTION_FINISH, removeRight);
			top.visible = false;
			tweenArray[4] = new Tween(right, 'y', Regular.easeOut, 0, sideSize, sideSpeed, true);
			tweenArray[5] = new Tween(right, 'height', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[5].addEventListener(TweenEvent.MOTION_FINISH, removeBottom);
		}
		
		private function removeBottom(e:TweenEvent):void{
			right.visible = false;
			tweenArray[6] = new Tween(bottom, 'width', Regular.easeOut, sideSize, 0, sideSpeed, true);
			tweenArray[6].addEventListener(TweenEvent.MOTION_FINISH, restartAnimation);
		}
		
		private function restartAnimation(e:TweenEvent):void{
			bottom.visible = false;
			tweenArray[6].removeEventListener(TweenEvent.MOTION_FINISH, restartAnimation);			
			this.startAnimation();
		}
		
		private function clearClips():void{
			left.graphics.clear();
			right.graphics.clear();
			top.graphics.clear();
			bottom.graphics.clear();
		}
			
	}
	
}