package projectbrowser{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	import flash.display.Loader;
	import flash.net.URLRequest;
	import flash.events.Event;
	import util.*;
	
	public class ProjectDetail extends MovieClip{
		
		private var defaultContainer:MovieClip;
		private var mouseOverContainer:MovieClip;
		private var projText:TextField;
		private var projTextFormat:TextFormat;
		private var verdana:Font = new VFont();
		private var nonDefaultText:String;
		private var defText:String;
		private var curMouseOver:Number = -1;
		private var curDefaultPath:String
		
		public function ProjectDetail():void{
			
			defaultContainer = new MovieClip();
			mouseOverContainer = new MovieClip();
			projTextFormat = new TextFormat();
			projTextFormat.size = 12;
			projTextFormat.font = verdana.fontName;
			projTextFormat.color = 0xCCCCCC;
			projText = new TextField();
			projText.defaultTextFormat = projTextFormat;
			projText.autoSize = TextFieldAutoSize.LEFT;
			projText.width = 255;
			projText.x = 10;
			this.addChild(defaultContainer);
			this.addChild(mouseOverContainer);
			this.addChild(projText);
			
		}
	
		public function setDefault(path:String, desc:String, projID:Number):void{
			//add the project details at this path and description to the default container
			//load the img at the path into the non-default container and show it
			//first, remove the old non-default content, if there was any
			curDefaultPath = path;
			if(defaultContainer.getChildByName('default')){
				defaultContainer.removeChild(defaultContainer.getChildByName('default'));
			}
			if(defaultContainer.getChildByName('defText')){
				defaultContainer.removeChild(defaultContainer.getChildByName('defText'));
			}
			if(projID == curMouseOver){
				if(mouseOverContainer.getChildByName('nonDefault')){
					var tmp:Loader = Loader(mouseOverContainer.getChildByName('nonDefault'));
					tmp.name = 'default';
					defaultContainer.addChild(tmp);
				}
				if(mouseOverContainer.getChildByName('nonDefText')){
					var tmp2:TextField = TextField(mouseOverContainer.getChildByName('nonDefText'));
					tmp2.name = 'defText';
					tmp2.width = 255;
					defaultContainer.addChild(tmp2);
				}
				mouseOverContainer.visible = false;
				defaultContainer.visible = true;
			}
			else{
				//set the non-default text var
				if(desc != "none"){
					defText = desc;
					var myFont:Font = new Font1();
					/* Create a new TextFormat object, and set the font property to the myFont
   					object's fontName property. */
					var myFormat:TextFormat = new TextFormat();
					myFormat.font = myFont.fontName;
					myFormat.size = 12;
	
					/* Create a new TextField object, assign the text format using the 
  					defaultTextFormat property, set the embedFonts property to true, and
 					 set the antiAliasType property to "advanced".*/
					var myTextField:TextField = new TextField();
					myTextField.width = 255;
					myTextField.autoSize = TextFieldAutoSize.LEFT;
					myTextField.defaultTextFormat = myFormat;
					myTextField.embedFonts = true;
					myTextField.antiAliasType = AntiAliasType.ADVANCED;
					myTextField.htmlText = desc;
					myTextField.name = 'defText';
					myTextField.x = 10;
					myTextField.wordWrap = true;
					myTextField.multiline = true;
					defaultContainer.addChild(myTextField);
				}
				mouseOverContainer.visible = false;
				defaultContainer.visible = true;
				if(path != "none"){
					if(defaultContainer.getChildByName('ldMovie')){
						defaultContainer.removeChild(defaultContainer.getChildByName('ldMovie'));
					}
					var ldMovie:LoadingMovie = new LoadingMovie(15, "Loading");
					ldMovie.x = (275 - ldMovie.width)/2;
					ldMovie.y = 20;
					ldMovie.name = 'ldMovie';
					ldMovie.startAnimation();
					defaultContainer.addChild(ldMovie);
					var loader:Loader = new Loader();
					loader.name = 'default';
					loader.contentLoaderInfo.addEventListener(Event.COMPLETE, finishDefaultUpLoad);
					loader.load(new URLRequest(MovieClip(root).browserImgBasePath + path));
					defaultContainer.addChild(loader);
				}
			}
			
		}
		
		private function finishDefaultUpLoad(e:Event):void{
			defaultContainer.removeChild(defaultContainer.getChildByName('ldMovie'));
			var ld:Loader = Loader(defaultContainer.getChildByName('default'));
			var txt:TextField = TextField(defaultContainer.getChildByName('defText'));
			txt.y = ld.y + ld.height + 5;
		}
		
		public function showDefault():void{
			//we're not moused over anything
			curMouseOver = -1;
			//make the default visible and hide the other
			mouseOverContainer.visible = false;
			defaultContainer.visible = true;
		}
		
		public function loadAndShow(path:String, desc:String, projID:Number):void{
			if(path != curDefaultPath){
				curMouseOver = projID;
				//load the img at the path into the non-default container and show it
				//first, remove the old non-default content, if there was any
				if(mouseOverContainer.getChildByName('nonDefault')){
					mouseOverContainer.removeChild(mouseOverContainer.getChildByName('nonDefault'));
				}
				if(mouseOverContainer.getChildByName('nonDefText')){
					mouseOverContainer.removeChild(mouseOverContainer.getChildByName('nonDefText'));
				}
				var myTextField:TextField = new TextField();
				//set the non-default text var
				if(desc != "none"){
					nonDefaultText = desc;
					/*var descFormat:TextFormat = new TextFormat;
					descFormat.font = verdana.fontName;
					descFormat.size = 12;
					descFormat.color = 0x000000;
					var descText:TextField = new TextField();
					descText.embedFonts = true;
					descText.setTextFormat(descFormat);
					descText.autoSize = TextFieldAutoSize.LEFT;
					descText.text = desc;
					descText.name = 'nonDefText';
					descText.width = 255;
					descText.x = 10;
					descText.wordWrap = true;
					descText.multiline = true;
					mouseOverContainer.addChild(descText);
					*/
					var myFont:Font = new Font1();
					/* Create a new TextFormat object, and set the font property to the myFont
   					object's fontName property. */
					var myFormat:TextFormat = new TextFormat();
					myFormat.font = myFont.fontName;
					myFormat.size = 12;

				/* Create a new TextField object, assign the text format using the 
  				defaultTextFormat property, set the embedFonts property to true, and
 				 set the antiAliasType property to "advanced". */
					myTextField.width = 255;
					myTextField.autoSize = TextFieldAutoSize.LEFT;
					myTextField.defaultTextFormat = myFormat;
					myTextField.embedFonts = true;
					myTextField.antiAliasType = AntiAliasType.ADVANCED;
					myTextField.htmlText = desc;
					myTextField.name = 'nonDefText';
					myTextField.x = 10;
					myTextField.wordWrap = true;
					myTextField.multiline = true;
					myTextField.visible = false;
					mouseOverContainer.addChild(myTextField);

				}
				defaultContainer.visible = false;
				mouseOverContainer.visible = true;
				if(path != "none"){
					if(mouseOverContainer.getChildByName('ldMovie')){
						mouseOverContainer.removeChild(mouseOverContainer.getChildByName('ldMovie'));
					}
					//create loading movie
					var ldMovie:LoadingMovie = new LoadingMovie(15, 'Loading');
					ldMovie.x = (275-ldMovie.width)/2;
					ldMovie.y = 20;
					ldMovie.name = 'ldMovie';
					ldMovie.startAnimation();
					mouseOverContainer.addChild(ldMovie);
					var loader:Loader = new Loader();
					loader.name = 'nonDefault';
					loader.contentLoaderInfo.addEventListener(Event.COMPLETE, finishUpLoad);
					loader.load(new URLRequest(MovieClip(root).browserImgBasePath + path));
					mouseOverContainer.addChild(loader);
				}
				else{
					myTextField.visible = true;
				}
			}
			
		}
		
		private function finishUpLoad(e:Event):void{
			if(mouseOverContainer.getChildByName('ldMovie')){
				mouseOverContainer.removeChild(mouseOverContainer.getChildByName('ldMovie'));
			}
			var ld:Loader = Loader(mouseOverContainer.getChildByName('nonDefault'));
			var txt:TextField = TextField(mouseOverContainer.getChildByName('nonDefText'));
			txt.visible = true;
			txt.y = ld.y + ld.height + 5;
			txt.width = 255;
		}
		
		
	}
	
	
}