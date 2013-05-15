package projectbrowser{
	
	import flash.display.MovieClip;
	import flash.display.Sprite;
	import flash.net.URLRequest;
	import flash.net.navigateToURL;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFieldType;
	import flash.text.Font;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import fl.containers.ScrollPane;
	import fl.containers.ScrollPaneNoBorder;
	import flash.events.Event;
	import flash.events.MouseEvent;
	import flash.events.ProgressEvent;
	import flash.net.Responder;
	import flash.net.NetConnection;
	import blast.BlastPanel;
	import detail.BottomPanel;
	import sequencepanel.ExpandEvent;
	import util.*;
	import flash.utils.Timer;
	
	public class CustomFasta extends MovieClip{
		
		private var fArea:CustomFastaArea;
		private var fBorder:Sprite;
		private var customSetBlast:CustomSetBlast;
		private var customDownload:CustomDownload;
		private var seqs:Array;
		private var currentSeqIndex:Number = 0;
		private var headerText:TextField;
		private var headerFormat:TextFormat = new TextFormat();
		private var vScrollBar:VScrollBar;
		private var seqPane:ScrollPaneNoBorder = new ScrollPaneNoBorder();
		private var seqContainer:MovieClip = new MovieClip();
		private var curY = 0;
		private var dBox:DownloadBox;
		private var downloadBar:DownloadBar;
		private var fileName:String = "";
		private var rollOverTimer:Timer = new Timer(100);
		
		public function CustomFasta(){
			//add the options for the Custom Fasta
			customDownload = new CustomDownload();
			customDownload.x = 20;
			customDownload.y = 0;
			customDownload.buttonMode = true;
			customDownload.name = 'download';
			customDownload.addEventListener(MouseEvent.ROLL_OVER, startTimer);
			customDownload.addEventListener(MouseEvent.CLICK, startDownload);
			this.addChild(customDownload);
			customSetBlast = new CustomSetBlast();
			customSetBlast.y = customDownload.y;
			customSetBlast.x = customDownload.x + customDownload.width;
			customSetBlast.buttonMode = true;
			customSetBlast.name = 'blast';
			customSetBlast.addEventListener(MouseEvent.ROLL_OVER, startTimer);
			customSetBlast.addEventListener(MouseEvent.CLICK, setCustomBlast);
			this.addChild(customSetBlast);
			seqContainer = new MovieClip();
			//make the sequence container
			seqPane.x = 10;
			seqPane.y = customSetBlast.y + customSetBlast.height + 10;
			seqPane.source = seqContainer;	
			seqPane.setSize(530, 180);
			seqPane.update();
			this.addChild(seqPane);
			seqs = new Array();
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
				if(mouseTarget.name == 'download'){
					aboutText.text = "download custom FASTA file";
				}
				else if(mouseTarget.name == 'blast'){
					aboutText.text = "set custom FASTA as BLAST query sequences";
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
				aboutBox.x = this.stage.mouseX;
				aboutBox.y = this.stage.mouseY + aboutBox.height;
				this.parent.parent.parent.addChild(aboutBox);
				mouseTarget.addEventListener(MouseEvent.MOUSE_OUT, destroyAbout);
			}
		}
		
		private function destroyAbout(e:MouseEvent):void{
			this.parent.parent.parent.removeChild(this.parent.parent.parent.getChildByName('aboutBox'));
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_OUT, destroyAbout);	
		}
		
		private function stopTimer(e:MouseEvent):void{
			e.currentTarget.removeEventListener(Event.ENTER_FRAME, checkTimer);
			e.currentTarget.removeEventListener(MouseEvent.MOUSE_OUT, stopTimer);
			rollOverTimer.reset();
		}
		
		public function getSeqs():Array{
			return seqs;			
		}
		
		
		private function startDownload(e:MouseEvent):void{
			if(seqs.length > 0){
				dBox = new DownloadBox();
				dBox.x = customSetBlast.x + customSetBlast.width + 5;;
				dBox.y = customSetBlast.y;
				this.addChild(dBox);
				dBox.gotoAndPlay(1);
				dBox.mouseChildren = false;
				dBox.buttonMode = true;
				var responder:Responder = new Responder(doDownload);
				var connection:NetConnection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				for(var i = 0; i < seqs.length; i++){
					trace(seqs[i][2]);
				}
				if(fileName == ""){
					connection.call("BLAST.downloadCustom", responder, seqs);
				}
				else{
					connection.call("BLAST.doNothing", responder);
				}
			}
		}
		
		private function doDownload(res:String):void{
			if(fileName == ""){
				fileName = res;
			}
			dBox.gotoAndStop(61);
			dBox.addEventListener(MouseEvent.CLICK, openLinkToFile);
			//create the file reference
			/*fastaRef = new FileReference();
			fastaRef.download(new URLRequest('http://10.41.128.72/seq_view/tmp/' + res), res);
			//create the download bar
			downloadBar = new DownloadBar(200, 15, 0xCCCCCC, 0x0099CC);
			downloadBar.y = fBorder.y + fBorder.height + 3;
			downloadBar.x = fBorder.x;
			this.addChild(downloadBar);
			//add the listeners
			fastaRef.addEventListener(ProgressEvent.PROGRESS, adjustDownloadBar);
			fastaRef.addEventListener(Event.COMPLETE, cleanUPDownload);
			fastaRef.addEventListener(Event.CANCEL, cleanUPDownload);
			*/
		}
		
		private function openLinkToFile(e:MouseEvent):void{
			dBox.removeEventListener(MouseEvent.CLICK, openLinkToFile);
			this.removeChild(dBox);
			navigateToURL(new URLRequest("../../../seq_view/tmp/" + fileName), "_blank");
		}
		
		private function cleanUPDownload(e:Event):void{
			/*//remove the downloadBar
			this.removeChild(downloadBar);
			//remove listeners
			fastaRef.removeEventListener(ProgressEvent.PROGRESS, adjustDownloadBar);
			fastaRef.removeEventListener(Event.COMPLETE, cleanUPDownload);
			fastaRef.removeEventListener(Event.CANCEL, cleanUPDownload);
			//remove the temporary file
			var responder = new Responder(removeComplete);
			var connection = new NetConnection();
			connection.connect(.gateway);
			connection.call("FlashDownload.removeTemp", responder, fileName);
			*/
		}
		
		private function adjustDownloadBar(e:ProgressEvent):void{

		}
		
		private function setCustomBlast(e:MouseEvent):void{
			//get the blast panel
			var bPanel:BlastPanel = BlastPanel(this.parent.parent.getChildByName('blast_panel'));
			//set the custom fasta as BLAST query
			bPanel.customCheck.selected = true;
			bPanel.fileCheck.selected = false;
		}
		
		private function removeComplete(res:String):void{
			trace(res);
		}
		
		public function addSequence(seqID, projectID, desc, defaultType){
			//invalidate the fileName
			fileName = "";
			//create a new entry in the array
			seqs[currentSeqIndex] = [seqID, projectID, defaultType];
			//display the newest sequence
			var tmp:CFastaSeq = new CFastaSeq(projectID, seqID, desc, seqContainer, currentSeqIndex, this);
			//add the listener to allow users to click this
			tmp.seqTextContainer.addEventListener(MouseEvent.CLICK, loadSeqDetail);
			//turn buttonmode on
			tmp.y = curY;
			//add the combobox
			if(defaultType == "NT"){
				createComboBox(currentSeqIndex, curY, 0);
			}
			else{
				createComboBox(currentSeqIndex, curY, 1);
			}
			seqContainer.addChild(tmp);
			seqPane.update();
			curY += tmp.height + 5;
			currentSeqIndex++;
			
		}
		
		private function loadSeqDetail(e:MouseEvent):void{
			//get the bottom panel
			var tmp:BottomPanel = BottomPanel(this.parent);
			//get the cFastaSeq
			var cFast:CFastaSeq = CFastaSeq(e.target.parent);
			//call sequence display method
			tmp.displaySequenceDetail(cFast.seqID, cFast.projID, cFast.desc, seqs[cFast.addID]);
		}
		
		private function createComboBox(index:Number, yVal:Number, select:Number):void{
			var comboBox:CustomComboBox = new CustomComboBox(index);
			comboBox.addItem({label:"NT", data:"NT"});
			comboBox.addItem({label:"AA", data:"AA"});
			comboBox.x = 25;
			comboBox.y = yVal; 
			comboBox.name = 'comboBox' + index;
			comboBox.selectedIndex = select;
			comboBox.width = 50;
			comboBox.addEventListener(Event.CHANGE, updateSeqType);
			seqContainer.addChild(comboBox);
		}
		
		private function updateSeqType(e:Event):void{
			var comboBox:CustomComboBox = CustomComboBox(e.target);
			if(comboBox.selectedIndex == 0){
				seqs[comboBox.customSeqID][2] = "NT";
			}
			else{
				seqs[comboBox.customSeqID][2] = "AA";
			}
		}
		
		public function removeSequence(seqIndex:Number):void{
			var tmp = seqContainer.getChildByName(String(seqIndex));
			seqContainer.removeChild(tmp);
			//remove the comboBox for this seq
			seqContainer.removeChild(seqContainer.getChildByName(String("comboBox" + seqIndex)));
			seqPane.update();
			//remove this element from the array
			restructureArray(seqIndex);
			this.dispatchEvent(new ExpandEvent(ExpandEvent.EXPAND_EVENT, tmp.addID, tmp.height +5));
			curY -= tmp.height + 5;
		}
		
		private function restructureArray(seqIndex:Number):void{
			for(var i = seqIndex; i < seqs.length - 1; i++){
				//move the seq array
				seqs[i] = seqs[i + 1];
			}
			seqs.pop();
			currentSeqIndex--;
			trace(seqs.length);
		}
		
		private function checkForScroll():void{
			if(seqContainer.height > 200){
				//resize and add the scroll bar
				//vScrollBar.resizeScroll(seqCotnainer.height);
				//this.addChild(VScrollBar);
			}
			else{
				//remove the scroll bar if it exists
				//if(this.contains(VScrollBar)){
					//this.removeChild(VScrollBar);
				//}
			}
		}
		
	}
	
	
}