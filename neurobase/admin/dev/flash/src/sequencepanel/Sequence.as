package sequencepanel{
	
	import flash.display.Sprite;
	import flash.display.MovieClip
	import flash.display.DisplayObject;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	import projectbrowser.*;
	import detail.*;
	import blast.BlastPanel;
	import util.*;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import fl.transitions.TweenEvent;
	import fl.transitions.Tween;
	import fl.transitions.easing.*;
	import flash.net.NetConnection;
	import flash.net.Responder;
	import flash.utils.Timer;

	
	public class Sequence extends Sprite{
		
		public var seqID;
		private var projectID;
		private var description;
		private var identifierText:TextField;
		private var rent:DisplayObject;
		private var expandAmount = 30;
		private var startY;
		private var endY;
		private var seqDetailWidth = 400;
		private var seqDetailHeight = 600;
		private var seqType;
		private var spacer = 5;
		private var rollOverTimer:Timer = new Timer(100);
		private var optionContainer:MovieClip;
		
		public function Sequence(ID:Number, proj:Number, desc:String, rent:DisplayObject, type:String){
			this.buttonMode = true;
			this.seqID = ID;
			this.projectID = proj;
			this.description = desc;
			//create the font
			var myFont:Font = new Font1();
			//create the fontformat for the sequence identifier
			var myFormat:TextFormat = new TextFormat();
			myFormat.font = myFont.fontName;
			myFormat.size = 12;
			this.identifierText = new TextField();
			identifierText.autoSize = TextFieldAutoSize.LEFT;
			identifierText.defaultTextFormat = myFormat;
			identifierText.embedFonts = true;
			identifierText.antiAliasType = AntiAliasType.ADVANCED;
			identifierText.width = seqDetailWidth;
			identifierText.wordWrap = true;
			identifierText.multiline = true;
			identifierText.text = desc;
			var textHolder:MovieClip = new MovieClip();
			textHolder.addChild(identifierText);
			textHolder.buttonMode = true;
			textHolder.mouseChildren = false;
			textHolder.addEventListener(MouseEvent.CLICK, startDetailLoad);
			this.addChild(textHolder);
			//rent.addEventListener(ExpandEvent.EXPAND_EVENT, onExpand);
			//this.addEventListener(MouseEvent.CLICK, dispatchExpand); 
			this.rent = rent;
			this.seqType = type;
			optionContainer = new MovieClip();
			//add option buttons
			var ntButton = new NucleotideButton();
			ntButton.x = 0;
			ntButton.addEventListener(MouseEvent.CLICK, displayNTSeq);
			ntButton.addEventListener(MouseEvent.MOUSE_OVER, startTimer);
			ntButton.name = 'ntButton';
			optionContainer.addChild(ntButton);
			var aaButton = new AminoAcidButton();
			aaButton.x = ntButton.x + ntButton.width + spacer;
			aaButton.addEventListener(MouseEvent.CLICK, displayAASeq);
			aaButton.addEventListener(MouseEvent.MOUSE_OVER, startTimer);
			aaButton.name = 'aaButton';
			optionContainer.addChild(aaButton);
			var blastButton = new BlastButton();
			blastButton.x = aaButton.x + aaButton.width + spacer;
			blastButton.addEventListener(MouseEvent.CLICK, addToCustomFasta);
			blastButton.addEventListener(MouseEvent.MOUSE_OVER, startTimer);
			blastButton.name = 'blastButton';
			optionContainer.addChild(blastButton);
			optionContainer.x = 535 - optionContainer.width - 5;
			this.addChild(optionContainer);
			//draw the background for the sequence
		}
		
		private function addToCustomFasta(e:MouseEvent):void{
			//check if there exists a custom Fasta
			var bPanel:BottomPanel = BottomPanel(this.rent.parent.getChildByName('bottom_panel'));
			bPanel.addCustomFastaSequence(seqID, projectID, description, seqType);
		}
		
		private function startDetailLoad(e:MouseEvent):void{
			if(this.seqType == 'AA'){
				displayAASeq(new MouseEvent(MouseEvent.CLICK));
			}
			else{
				displayNTSeq(new MouseEvent(MouseEvent.CLICK));
			}
		}
		
		private function startTimer(e:MouseEvent):void{
			rollOverTimer.start();
			e.currentTarget.addEventListener(Event.ENTER_FRAME, checkTimer);
			e.currentTarget.addEventListener(MouseEvent.MOUSE_OUT, stopTimer);
		}
		
		private function checkTimer(e:Event):void{
			if(rollOverTimer.currentCount > 4){
				//stop the timer
				rollOverTimer.reset();
				//remove the current mouse out listener
				e.currentTarget.removeEventListener(MouseEvent.MOUSE_OUT, stopTimer);
				e.currentTarget.removeEventListener(Event.ENTER_FRAME, checkTimer);
				//display the text field about the button
				var mouseTarget:MovieClip = MovieClip(e.currentTarget);
				var aboutText:TextField = new TextField();
				aboutText.autoSize = TextFieldAutoSize.LEFT;
				if(mouseTarget.name == 'ntButton'){
					aboutText.text = "view nucleotide sequence";
				}
				else if(mouseTarget.name == 'aaButton'){
					aboutText.text = "view amino acid sequence";
				}
				else if(mouseTarget.name == 'blastButton'){
					aboutText.text = "add to Custom Fasta File";
				}
				//draw the containing box
				var aboutBox:MovieClip = new MovieClip();
				var aboutBoxWidth:Number = aboutText.width + 10;
				var aboutBoxHeight:Number = aboutText.height + 5;
				aboutBox.graphics.beginFill(0xFEFCCA, 1);
				aboutBox.graphics.drawRect(0, 0, aboutBoxWidth, aboutBoxHeight);
				aboutBox.graphics.endFill();
				aboutBox.graphics.beginFill(0x000000, 1);
				aboutBox.graphics.drawRect(0, 0, aboutBoxWidth, 1);
				aboutBox.graphics.drawRect(0, 0, 1, aboutBoxHeight);
				aboutBox.graphics.drawRect(0, aboutBoxHeight, aboutBoxWidth, 1);
				aboutBox.graphics.drawRect(aboutBoxWidth, 0, 1, aboutBoxHeight);
				aboutBox.graphics.endFill();
				aboutText.x = 5;
				aboutText.y = 2.5;
				aboutBox.addChild(aboutText);
				aboutBox.name ='aboutBox';
				aboutBox.x = this.stage.mouseX - aboutBox.width;
				aboutBox.y = this.stage.mouseY - aboutBox.height;
				this.rent.parent.parent.addChild(aboutBox);
				mouseTarget.addEventListener(MouseEvent.MOUSE_OUT, destroyAbout);
			}
		}
		
		private function destroyAbout(e:MouseEvent):void{
			this.rent.parent.parent.removeChild(this.rent.parent.parent.getChildByName('aboutBox'));
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_OUT, destroyAbout);	
		}
		
		private function stopTimer(e:MouseEvent):void{
			e.currentTarget.removeEventListener(Event.ENTER_FRAME, checkTimer);
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_OUT, stopTimer);
			rollOverTimer.reset();
		}
		
		private function displayNTSeq(e:MouseEvent):void{
			var dPanel:BottomPanel = BottomPanel(this.rent.parent.getChildByName('bottom_panel'));
			dPanel.displaySequenceDetail(this.seqID, this.projectID, this.description, "NT");
		}
		
		private function displayAASeq(e:MouseEvent):void{
			var dPanel:BottomPanel = BottomPanel(this.rent.parent.getChildByName('bottom_panel'));
			dPanel.displaySequenceDetail(this.seqID, this.projectID, this.description, "AA");
		}
		
		public function initYs(){
			startY = this.y;
			endY = this.y + expandAmount;
		}
		
				
		internal function displayDetails(e:MouseEvent){
			var bottomPanel:BottomPanel = BottomPanel(this.rent.parent.getChildByName('bottom_panel'));
			bottomPanel.displaySequenceDetail(this.seqID, this.projectID, this.description, this.seqType);
		}
		
		internal function querySeq(e:MouseEvent){
			var responder:Responder = new Responder(showSeq);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query = "SELECT seq_id, description, sequence, type FROM " + this.projectID + "_sequences WHERE seq_id = '" + this.seqID + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function showSeq(res:Array){
			var totalHeight:Number = 0;
			var space:Number = 30;
			var sequence:String = res[0]['sequence'];
			var type = res[0]['type'];
			var seqText:TextField = new TextField();
			var seqBackground:MovieClip = new MovieClip();
			seqBackground.name = 'seqBackground';
			parent.addChild(seqBackground);			
			/*verticalScroll.direction = "vertical";
			verticalScroll.setSize(seqText.width, seqText.height);
			verticalScroll.move(seqText.x + seqText.width, seqText.y);
			verticalScroll.scrollTarget = seqText;
			*/
			var closeText:TextClip = new TextClip('x close');
			closeText.x = seqDetailWidth - 5 - closeText.width;
			closeText.y = 0;
			closeText.addEventListener(MouseEvent.CLICK, closeSeqDetail);
			seqText.text = ">" + res[0]['description'] + "\n\n" + sequence;
			seqText.x = 10;
			seqText.y = closeText.y + closeText.height + space;
			seqText.width = seqDetailWidth - 40;
			seqText.multiline = true;
			seqText.wordWrap = true;
			seqText.autoSize = TextFieldAutoSize.LEFT;
			totalHeight = closeText.height + space + seqText.height + space;
			seqBackground.addChild(closeText);
			seqBackground.addChild(seqText);
			seqBackground.graphics.beginFill(0xFFFFFF, 1);
			seqBackground.graphics.drawRect(0, 0, seqDetailWidth, totalHeight);
			seqBackground.graphics.endFill();
			seqBackground.name = 'seqBackground'
			seqBackground.graphics.beginFill(0x000000, 1);
			seqBackground.graphics.drawRect(0, 0, 1, totalHeight);
			seqBackground.graphics.drawRect(seqDetailWidth - 1, 0, 1, totalHeight);
			seqBackground.graphics.drawRect(0, 0, seqDetailWidth - 1, 1);
			seqBackground.graphics.drawRect(0, totalHeight, seqDetailWidth - 1, 1);
			seqBackground.graphics.endFill();
			seqBackground.y = this.y;
			SequencePanel(parent.parent.parent).resizeScroll('expand', 0);
			//check to see if there the other version of this sequence exists in the DB
			var responder:Responder = new Responder(checkOtherType);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query;
			if(type == "AA"){
				//check for a nucleotide sequence
				query = "SELECT seq_id FROM " + this.projectID + "_sequences WHERE description = '" + this.description + "' and type='NT'";
			
			}
			else if(type == "NT"){
				//check for an amino acid sequence
				query = "SELECT seq_id FROM " + this.projectID + "_sequences WHERE description = '" + this.description + "' and type='AA'";
			
			}
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function checkOtherType(res:Array){
			if(res.length > 0){
				//display option to view the other type of sequence
			}
		}		
		private function closeSeqDetail(e:MouseEvent){
			var tmp = MovieClip(parent.getChildByName('seqBackground'));
			parent.removeChild(tmp);
			SequencePanel(parent.parent.parent).resizeScroll('contract', tmp.height);
		}
		
		private function dispatchExpand(e:Event){
			this.removeEventListener(MouseEvent.CLICK, dispatchExpand);
			this.addEventListener(MouseEvent.CLICK, dispatchContract);
			this.parent.dispatchEvent(new ExpandEvent('expandevent', seqID, expandAmount));			
		}
		
		private function dispatchContract(e:Event){
	
			this.removeEventListener(MouseEvent.CLICK, dispatchContract);
			this.addEventListener(MouseEvent.CLICK, dispatchExpand);
			this.parent.dispatchEvent(new ExpandEvent('expandevent', seqID, (expandAmount * -1)));
		}
		
		internal function addToBlast(e:MouseEvent){
			var responder:Responder = new Responder(addQuerySeq);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query = "SELECT seq_id, description, sequence FROM " + this.projectID + "_sequences WHERE seq_id = '" + this.seqID + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		internal function addQuerySeq(res:Array){
			BlastPanel(this.parent.parent.parent.parent.getChildByName('blast_panel')).addQuerySeq(res[0]['description'], res[0]['sequence']);
		}
		
	}
	
	
}