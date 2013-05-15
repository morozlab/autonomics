package blast{
	
	import flash.display.MovieClip;
	import flash.display.DisplayObject;
	import fl.containers.ScrollPane;
	import flash.events.MouseEvent;
	import flash.events.Event;
	import flash.display.LoaderInfo;
	import fl.transitions.Tween;
	import fl.transitions.easing.*;
	import fl.transitions.TweenEvent;
	import flash.display.Loader;
	import flash.net.URLRequest;
	import flash.net.URLLoader;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.text.Font;
	import flash.net.NetConnection;
	import flash.net.Responder;
	import util.*;
	
	
	public class BlastDisplay extends MovieClip{
		
		private var displayWidth:Number;
		private var blackBG:MovieClip = new MovieClip();
		private var alphaTween:Tween;
		private var utilTween:Tween;
		private var whitePart:MovieClip = new MovieClip();
		private var brButton:MovieClip = new MovieClip();
		private var quantButton:MovieClip = new MovieClip();
		private var graphicsButton:MovieClip = new MovieClip();
		private var closeB:MovieClip = new MovieClip();
		private var blastContainer:MovieClip = new MovieClip();
		private var quantContainer:MovieClip = new MovieClip();
		private var graphicsContainer:MovieClip = new MovieClip();
		private var contentContainer:MovieClip = new MovieClip();
		private var contentMask:MovieClip = new MovieClip();
		private var graphicsArray:Array = new Array();
		private var blastText:Array = new Array();
		private var blastID:Number;
		private var loadedVars:Array;
		private var curGraphicIndex = 0;
		private var curGraphicY = 0;
		
		public function BlastDisplay(w:Number){
			this.displayWidth = w;
			this.name ='blast_display';
			
		}
		
		public function initDisplay(blast:Number):void{
			//add the display elements
			initDisplayElements();
			blackBG.name = 'blackBG';
			blackBG.graphics.beginFill(0x000000, 1);
			blackBG.graphics.drawRect(0, 0, 1180, 1200);
			blackBG.graphics.endFill();
			whitePart.graphics.beginFill(0xFFFFFF, 1);
			whitePart.graphics.drawRect(0, 0, displayWidth, 900);
			whitePart.graphics.endFill();
			whitePart.x = (1180-800)/2;
			whitePart.name = 'whitePart';
			this.addChild(blackBG);
			//alpha tween on the black BG
			alphaTween = new Tween(blackBG, 'alpha', Regular.easeOut, 0, 1, .5, true);
			contentContainer.x = (1180 - 800)/2;
			contentContainer.name = 'contentContainer';
			contentMask.graphics.beginFill(0, 0);
			contentMask.graphics.drawRect(0, 0, 900, 900);
			contentMask.graphics.endFill();
			contentMask.x = contentContainer.x;
			contentContainer.mask = contentMask;
			this.addChild(contentContainer);
			blastID = blast;
			alphaTween.addEventListener(TweenEvent.MOTION_FINISH, drawWhite);
		}
		
		private function initDisplayElements():void{
			blackBG = new MovieClip();
			whitePart = new MovieClip();
			contentContainer = new MovieClip();
			brButton = new MovieClip();
			quantButton = new MovieClip();
			graphicsButton = new MovieClip();
			closeB = new MovieClip();
			blastContainer = new MovieClip();
			quantContainer = new MovieClip();
			graphicsContainer = new MovieClip();
			reinitArray();
		}
		
		private function reinitArray():void{
			loadedVars = new Array(0, 0, 0);
		}
		
		private function drawWhite(e:TweenEvent):void{
			blackBG.addChild(whitePart);
			utilTween = new Tween(whitePart, 'width', Regular.easeOut, 0, displayWidth, .5, true);
			utilTween.addEventListener(TweenEvent.MOTION_FINISH, addNavigation);
		}
		
		private function addNavigation(e:TweenEvent):void{
			var brText:BlastResultButton = new BlastResultButton();
			brButton.x = whitePart.x - brText.width - 5;
			brButton.y = 200;
			brButton.buttonMode = true;
			brButton.addEventListener(MouseEvent.ROLL_OVER, brMouseOver);
			brButton.addEventListener(MouseEvent.CLICK, handleBlastClick);
			brText.name = 'out';
			brButton.addChild(brText);
			this.addChild(brButton);
			var quantText:QuantButton = new QuantButton();
			quantButton.x = whitePart.x - brText.width - 5;
			quantButton.y = brButton.y + brButton.height + 3;
			quantButton.buttonMode = true;
			quantButton.addEventListener(MouseEvent.CLICK, handleQuantClick);
			quantButton.addEventListener(MouseEvent.ROLL_OVER, quantMouseOver);
			quantText.name = 'out';
			quantButton.addChild(quantText);
			this.addChild(quantButton);
			var graphicsText:GraphicsButton = new GraphicsButton();
			graphicsButton.x = whitePart.x - brText.width - 5;
			graphicsButton.y = quantButton.y + quantButton.height + 3;
			graphicsButton.buttonMode = true;
			graphicsButton.addEventListener(MouseEvent.ROLL_OVER, graphicMouseOver);
			graphicsButton.addEventListener(MouseEvent.CLICK, handleGraphicsClick);
			graphicsText.name = 'out';
			graphicsButton.addChild(graphicsText);
			var closeText:BlastClose = new BlastClose();
			closeB.x = whitePart.x - closeText.width - 5;
			closeB.y = graphicsButton.y + graphicsButton.height + 3;
			closeB.buttonMode = true;
			closeB.addEventListener(MouseEvent.ROLL_OVER, closeMouseOver);
			closeB.addEventListener(MouseEvent.CLICK, doClose);
			closeText.name = 'out';
			closeB.addChild(closeText);
			this.addChild(closeB);
			this.addChild(graphicsButton);
			this.loadBlastResult();
			
		}
		
		private function closeMouseOver(e:MouseEvent):void{
			closeB.removeEventListener(MouseEvent.ROLL_OVER, closeMouseOver);
			while(closeB.getChildByName('out')){
				closeB.removeChild(closeB.getChildByName('out'));
			}
			var closeBOver:BlastCloseOver = new BlastCloseOver();
			closeBOver.name = 'over';
			closeB.addChild(closeBOver);
			closeB.addEventListener(MouseEvent.ROLL_OUT, closeMouseOut);
		}
		
		private function closeMouseOut(e:MouseEvent):void{
			closeB.removeEventListener(MouseEvent.ROLL_OVER, closeMouseOut);
			while(closeB.getChildByName('over')){
				closeB.removeChild(closeB.getChildByName('over'));
			}
			var closeBOver:BlastClose = new BlastClose();
			closeBOver.name = 'out';
			closeB.addChild(closeBOver);
			closeB.addEventListener(MouseEvent.ROLL_OVER, closeMouseOver);
		}
		
		private function doClose(e:MouseEvent):void{
			this.clearDisplay();
			//start tween of white thing
			utilTween = new Tween(whitePart, 'width', Regular.easeOut, whitePart.width, 0, .5, true);
			utilTween.addEventListener(TweenEvent.MOTION_FINISH, exitAlphaTween);
		}
		
		private function exitAlphaTween(e:TweenEvent):void{
			utilTween.removeEventListener(TweenEvent.MOTION_FINISH, exitAlphaTween);
			utilTween = new Tween(blackBG, 'alpha', Regular.easeOut, 1, 0, .5, true);
			utilTween.addEventListener(TweenEvent.MOTION_FINISH, selfRemoval);
		}
		
		private function selfRemoval(e:TweenEvent):void{
			this.removeChild(brButton);
			this.removeChild(quantButton);
			this.removeChild(closeB);
			this.removeChild(graphicsButton);
			blackBG.removeChild(whitePart);
			this.removeChild(contentContainer);
			if(this.contains(blackBG)){
				this.removeChild(blackBG);
			}
			this.parent.removeChild(this);
		}
				
		private function graphicMouseOver(e:MouseEvent):void{
			graphicsButton.removeEventListener(MouseEvent.ROLL_OVER, graphicMouseOver);
			while(graphicsButton.getChildByName('out')){
				graphicsButton.removeChild(graphicsButton.getChildByName('out'));
			}
			var brOverText:GraphicsOverButton = new GraphicsOverButton();
			brOverText.name = 'over';
			graphicsButton.addChild(brOverText);
			graphicsButton.addEventListener(MouseEvent.ROLL_OUT, graphicMouseOut);
			
		}
		
		private function graphicMouseOut(e:MouseEvent):void{
			graphicsButton.removeEventListener(MouseEvent.ROLL_OUT, graphicMouseOut);
			while(graphicsButton.getChildByName('over')){
				graphicsButton.removeChild(graphicsButton.getChildByName('over'));
			}
			var brOverText:GraphicsButton = new GraphicsButton();
			brOverText.name = 'out';
			graphicsButton.addChild(brOverText);
			graphicsButton.addEventListener(MouseEvent.ROLL_OVER, graphicMouseOver);
		}
		
		private function brMouseOver(e:MouseEvent):void{
			brButton.removeEventListener(MouseEvent.ROLL_OVER, brMouseOver);
			while(brButton.getChildByName('out')){
				brButton.removeChild(brButton.getChildByName('out'));
			}
			var brOverText:BlastResultOver = new BlastResultOver();
			brOverText.name = 'over';
			brButton.addChild(brOverText);
			brButton.addEventListener(MouseEvent.ROLL_OUT, brMouseOut);
		}
		
		private function brMouseOut(e:MouseEvent):void{
			brButton.removeEventListener(MouseEvent.ROLL_OUT, brMouseOut);
			while(brButton.getChildByName('over')){
				brButton.removeChild(brButton.getChildByName('over'));
			}
			var brOverText:BlastResultButton = new BlastResultButton();
			brOverText.name = 'out';
			brButton.addChild(brOverText);
			brButton.addEventListener(MouseEvent.ROLL_OVER, brMouseOver);
		}
		
		private function quantMouseOver(e:MouseEvent):void{
			quantButton.removeEventListener(MouseEvent.ROLL_OVER, quantMouseOver);
			while(quantButton.getChildByName('out')){
				quantButton.removeChild(quantButton.getChildByName('out'));
			}
			var brOverText:QuantOverButton = new QuantOverButton();
			brOverText.name = 'over';
			quantButton.addChild(brOverText);
			quantButton.addEventListener(MouseEvent.ROLL_OUT, quantMouseOut);
		}
		
		private function quantMouseOut(e:MouseEvent):void{
			quantButton.removeEventListener(MouseEvent.ROLL_OUT, quantMouseOut);
			while(quantButton.getChildByName('over')){
				quantButton.removeChild(quantButton.getChildByName('over'));
			}
			var brOverText:QuantButton = new QuantButton();
			brOverText.name = 'out';
			quantButton.addChild(brOverText);
			quantButton.addEventListener(MouseEvent.ROLL_OVER, quantMouseOver);			
		}
//***************************************************** BLAST FUNCTIONS *********************************************//
		private function handleBlastClick(e:MouseEvent):void{
			this.clearDisplay();
			if(loadedVars[0] == 0){
				this.loadBlastResult();	
			}
			else{
				this.reloadBlastContainer();
			}
		}
		
		private function loadBlastResult():void{
			//load the text from the BLAST file
			var textLoader:URLLoader = new URLLoader();
			textLoader.addEventListener(Event.COMPLETE, parseBlastOutput);
			textLoader.load(new URLRequest("../../../seq_view/results/" + blastID + "/BLAST_Results/blastOutput.out"));
		}
		
		private function handleBlastText(e:Event):void{
			
		}
		
		private function parseBlastOutput(e:Event):void{
			var cFont:Courier = new Courier();
			var blastFormat:TextFormat = new TextFormat();
			blastFormat.font = cFont.fontName;
			blastFormat.size = 14;
			var textBox:TextField = new TextField();
			textBox.autoSize = TextFieldAutoSize.LEFT;
			textBox.width = 760;
			textBox.height = 700;
			textBox.x = 20;
			textBox.y = 10;
			textBox.name = 'blastText';
			textBox.defaultTextFormat = blastFormat;
			textBox.setTextFormat(blastFormat);
			textBox.text = e.target.data;
			blastContainer.addChild(textBox);
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = blastContainer;
			blastPane.name = 'blastPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
			//for(var i = 0; i < res.length; i++){
			//	textBox.appendText(res[i] + "\n");
			//}
			loadedVars[0] = 1;			
		}
	
		private function reloadBlastContainer():void{
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = blastContainer;
			blastPane.name = 'blastPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
		}

//****************************************** QUANTIFICATION FUNCTIONS ********************************************//

		private function handleQuantClick(e:MouseEvent):void{
			this.clearDisplay()
			if(loadedVars[1] == 0){
				this.loadQuantification();
			}
			else{
				this.reloadQuantContainer();
			}
		}
		
		private function loadQuantification():void{
			trace('moo');
			var responder = new Responder(doQuantLoad);
			var connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			connection.call("BLAST.loadQuantification", responder, blastID);
		}
		
		private function doQuantLoad(res:Array):void{
			trace('moo2');
			var vFont:Font = new Font1();
			var quantFormat:TextFormat = new TextFormat();
			quantFormat.font = vFont.fontName;
			quantFormat.size = 12;
			//add the quantification container
			this.contentContainer.addChild(quantContainer);
			var curY = 0;
			//print rows
			for(var i =0; i < res.length; i++){
				var line:String = res[i];
				var splitLine:Array = line.split(/\t/);
				//print columns
				var curX = 20;
				var maxHeight:Number = 0;
				for(var j = 0; j < splitLine.length; j++){
					var quantText:TextField = new TextField();
					quantText.defaultTextFormat = quantFormat;
					quantText.embedFonts = true;
					quantText.wordWrap = true;
					quantText.multiline = true;
					quantText.x = curX;
					quantText.y = curY;
					quantText.text = splitLine[j];
					
					if(j == 0){
						quantText.autoSize = TextFieldAutoSize.LEFT;
						quantText.width = 345;
						curX += 345;
						maxHeight = quantText.height;
					}
					else{
						quantText.autoSize = TextFieldAutoSize.RIGHT;
						quantText.width = 145;
						curX += 145;
					}
					quantContainer.addChild(quantText);
					if(j == (splitLine.length - 1)){
						if(i == 1){
							curY += 50;
						}
						else{
							curY += maxHeight + 5;
						}
					}
				}
			}
			loadedVars[1] = 1;
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = quantContainer;
			blastPane.name = 'quantPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
		}
		
		private function reloadQuantContainer():void{
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = quantContainer;
			blastPane.name = 'quantPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
		}

//******************************************* GRAPHICS FUNCTIONS *************************************************//

	private function handleGraphicsClick(e:MouseEvent):void{
		this.clearDisplay();
		if(loadedVars[2] == 0){
			curGraphicIndex = 0;
			curGraphicY = 0;
			this.getGraphicFilenames();
		}
		else{
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = graphicsContainer;
			blastPane.name = 'graphicsPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
		}
	}
	
	private function getGraphicFilenames():void{
		var responder:Responder = new Responder(displayGraphics);
		var connection:NetConnection = new NetConnection();
		connection.connect(MovieClip(root).gateway);
		connection.call("BLAST.getGraphicFilenames", responder, blastID);
		
	}
	
	private function displayGraphics(res:Array):void{
		if(res.length > 0){//save the graphics array
			this.graphicsArray = res;
			//add the graphics container
			this.contentContainer.addChild(graphicsContainer);
			//load the first image
			this.loadNextGraphic(curGraphicIndex);
		}
	}
	
	private function loadNextGraphic(index:Number):void{
		trace(graphicsArray[index]);
		var loader:Loader = new Loader();
		var req:URLRequest = new URLRequest(MovieClip(root).basePath + '/seq_view/results/' + blastID + '/Graphics/' + graphicsArray[curGraphicIndex]);
		loader.x = 50;
		loader.y = curGraphicY;
		loader.contentLoaderInfo.addEventListener(Event.COMPLETE, handleGraphicLoad);
		loader.load(req);
		graphicsContainer.addChild(loader);
		curGraphicIndex = index + 1;
	}
	
	private function handleGraphicLoad(e:Event):void{
		if(curGraphicIndex < graphicsArray.length){
			curGraphicY += LoaderInfo(e.currentTarget).height + 20;
			loadNextGraphic(curGraphicIndex);
		}		
		else{
			this.loadedVars[2] = 1;
			var blastPane:ScrollPane = new ScrollPane();
			blastPane.source = graphicsContainer;
			blastPane.name = 'graphicsPane';
			blastPane.setSize(800, 900);
			blastPane.verticalScrollPolicy = "auto";
			blastPane.update();
			this.contentContainer.addChild(blastPane);
		}
	}
	
	

//******************************************* OTHER STUFF ********************************************************//
		
		private function clearDisplay():void{
			
			if(this.contentContainer.getChildByName('blastPane')){
				this.contentContainer.removeChild(this.contentContainer.getChildByName('blastPane'));
			}
			if(this.contentContainer.getChildByName('quantPane')){
				this.contentContainer.removeChild(this.contentContainer.getChildByName('quantPane'));
			}
			if(this.contentContainer.getChildByName('graphicsPane')){
				this.contentContainer.removeChild(this.contentContainer.getChildByName('graphicsPane'));
			}
			if(this.getChildByName('scrollBar')){
				this.removeChild(this.getChildByName('scrollBar'));
			}
			
		}
		
		private function resetDisplay():void{
			blastContainer = new MovieClip();
			quantContainer = new MovieClip();
			graphicsContainer = new MovieClip();
			loadedVars = new Array(0, 0, 0);
			graphicsArray = new Array();
			curGraphicIndex = 0;
			curGraphicY = 0;
		}
		
		private function resetBG():void{
			brButton = new MovieClip();
			quantButton = new MovieClip();
			graphicsButton = new MovieClip();
			closeB = new MovieClip();
			if(this.getChildByName('blackBG')){
				this.removeChild(blackBG);
			}
			blackBG.graphics.clear();
			if(blackBG.getChildByName('whitePart')){
				blackBG.removeChild(whitePart);
			}
			whitePart.graphics.clear();
		}
		
		private function loadScrollBar(h:Number, target:DisplayObject):void{
			if(h >= contentMask.height){
				var blastScroll:VScrollBar = new VScrollBar(target, h, 900);
				blastScroll.x = contentContainer.x + 800;
				blastScroll.y = 0;
				blastScroll.name = 'scrollBar';
				this.addChild(blastScroll);
			}
		}
	}
	
	
}