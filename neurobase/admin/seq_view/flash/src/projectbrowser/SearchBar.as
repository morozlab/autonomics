package projectbrowser{
	
	import flash.display.Sprite;
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFieldType;
	import flash.events.MouseEvent;
	import flash.events.KeyboardEvent;
	import fl.transitions.Tween;
	import fl.transitions.easing.*;
	import fl.transitions.TweenEvent;
		
	public class SearchBar extends Sprite{
		
		public var titleText = new TextField();
		public var searchText:TextField = new TextField();
		private var top:MovieClip = new MovieClip();
		private var bottom:MovieClip = new MovieClip();
		private var left:MovieClip = new MovieClip();
		private var right:MovieClip = new MovieClip();
		var tTween:Tween;
		var bTween:Tween;
		var bxTween:Tween;
		var lTween:Tween;
		var rTween:Tween;
		var rightYTween:Tween;
		//public var clickBox = new SearchIcon();
		
		public function SearchBar():void{
			//make four movie clips with black background
			top.graphics.beginFill(0x000000, 1);
			top.graphics.drawRect(0, 0, 1, 1);
			top.graphics.endFill();
			bottom.graphics.beginFill(0x000000, 1);
			bottom.graphics.drawRect(0, 0, 1, 1);
			bottom.graphics.endFill();
			bottom.y = 20;
			bottom.x = 100;
			left.graphics.beginFill(0x000000, 1);
			left.graphics.drawRect(0, 0, 1, 1);
			left.graphics.endFill();
			right.graphics.beginFill(0x000000, 1);
			right.graphics.drawRect(0, 0, 1, 1);
			right.graphics.endFill();
			right.x = 100;
			right.y = 20;
			//searchText.border = true;
			//searchText.borderColor = 0x000000;
			searchText.width = 100;
			searchText.type = TextFieldType.INPUT;
			searchText.height = 19;
			searchText.addEventListener(KeyboardEvent.KEY_DOWN, keyExecSearch);
			searchText.text = "";
			//clickBox.x = searchText.x + searchText.width + 3;
			//clickBox.y = (searchText.height - clickBox.height + 2) /2
			//clickBox.addEventListener(MouseEvent.CLICK, executeSearch);
			//this.addChild(clickBox);
			this.addChild(searchText);
		}
		
		public function startTween():void{
			this.addChild(top);
			tTween = new Tween(top, 'width', Regular.easeOut, 0, 100, .5, true);
			this.addChild(bottom);
			bTween = new Tween(bottom, 'width', Regular.easeOut, 0, 100, .5, true);
			bxTween = new Tween(bottom, 'x', Regular.easeOut, 100, 0, .5, true);
			this.addChild(left);
			lTween = new Tween(left, 'height', Regular.easeOut, 0, 20, .5, true);
			this.addChild(right);
			rTween = new Tween(right, 'height', Regular.easeOut, 0, 20, .5, true);
			rightYTween = new Tween(right, 'y', Regular.easeOut, 20, 0, .5, true);
			rightYTween.addEventListener(TweenEvent.MOTION_FINISH, enableClick);
		}
		
		private function enableClick(e:TweenEvent):void{
			ProjectDisplay(parent).enableSearchClick();
			rightYTween.removeEventListener(TweenEvent.MOTION_FINISH, enableClick);
		}
		
		private function executeSearch(e:MouseEvent):void{
			var tmp:ProjectDisplay = ProjectDisplay(this.parent);
			tmp.execSearch(searchText.text);
		}
		
		private function keyExecSearch(e:KeyboardEvent):void{
			if(e.charCode == 13){
				this.hideBoxAndSearch();
			}
		}
		
		internal function hideBoxAndSearch():void{
			searchText.visible = false;
			tTween = new Tween(top, 'width', Regular.easeOut, 100, 0, .5, true);
			bTween = new Tween(bottom, 'width', Regular.easeOut, 100, 0, .5, true);
			bxTween = new Tween(bottom, 'x', Regular.easeOut, 0, 100, .5, true);
			lTween = new Tween(left, 'height', Regular.easeOut, 20, 0, .5, true);
			rTween = new Tween(right, 'height', Regular.easeOut, 20, 0, .5, true);
			rightYTween = new Tween(right, 'y', Regular.easeOut, 0, 20, .5, true);
			//execute the search
			while(rightYTween.hasEventListener(TweenEvent.MOTION_FINISH)){
				rightYTween.removeEventListener(TweenEvent.MOTION_FINISH, hideSides);
			}
			rightYTween.addEventListener(TweenEvent.MOTION_FINISH, hideSides);
		}
		
		private function hideSides(e:TweenEvent):void{
			rightYTween.removeEventListener(TweenEvent.MOTION_FINISH, hideSides);
			this.removeChild(top);
			this.removeChild(bottom);
			this.removeChild(left);
			this.removeChild(right);
			ProjectDisplay(parent).disableSearchClick();
			if(searchText.text != ""){
				var tmp:ProjectDisplay = ProjectDisplay(this.parent);
				tmp.execSearch(searchText.text);
			}				
		}
		
		
		
	}
	
}