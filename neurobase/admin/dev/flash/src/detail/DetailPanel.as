package detail{
	
	import flash.display.MovieClip;
	import util.*;
	import flash.net.Responder;
	import flash.net.NetConnection;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	import flash.events.MouseEvent;
	import fl.containers.ScrollPane;
	import fl.containers.ScrollPaneNoBorder;
	import sequencepanel.*;
	import blast.BlastPanel;
	
	public class DetailPanel extends MovieClip{
		
		private var detailPanelWidth = 530;
		private var detailPanelHeight = 200;
		private var seqDetailPane:ScrollPaneNoBorder; 
		private var detailContainer:MovieClip = new MovieClip();
		private var seqID:Number;
		private var projectID:Number;
		private var desc:String;
		private var detailMask:MovieClip = new MovieClip();
		private var seqType:String;
		private var initContainerY:Number = 30;
		
		public function DetailPanel(){
			detailContainer.x = 0;
			detailContainer.y = 0;
			detailContainer.graphics.beginFill(0, 0);
			detailContainer.graphics.drawRect(0, 0, detailPanelWidth, detailPanelHeight);
			detailContainer.graphics.endFill();	
			this.name = 'detail_panel';
			seqDetailPane = new ScrollPaneNoBorder();
			seqDetailPane.source = detailContainer;
			seqDetailPane.setSize(detailPanelWidth + 20, detailPanelHeight);
			this.addChild(seqDetailPane);
		}
		
		public function displaySequenceDetail(seqID:Number, projectID:Number, description:String, type:String){
			this.seqID = seqID;
			this.projectID = projectID;
			this.desc = description;
			this.seqType = type;
			var responder:Responder = new Responder(showSeq);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query;
			if(seqType == "AA"){
				query = "SELECT seq_id, sb_id, description, AA_sequence FROM " + this.projectID + "_sequences WHERE seq_id= '" + seqID + "'";
			}
			else{
				query = "SELECT seq_id, sb_id, description, NT_sequence FROM " + this.projectID + "_sequences WHERE seq_id= '" + seqID + "'";
			}
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function showSeq(res:Array){
			//clear all of the things being displayed in the detail panel
			if(detailContainer.getChildByName('seqText')){
				detailContainer.removeChild(detailContainer.getChildByName('seqText'));
				if(this.getChildByName('changeText')){
					this.removeChild(this.getChildByName('changeText'));
				}
				if(this.getChildByName('blastClip')){
					this.removeChild(this.getChildByName('blastClip'));
				}
			}
			if(res.length > 0){
				var totalHeight:Number = 0;
				var space:Number = 0;
			var mod:String;
			if(seqType == "AA"){
				mod = "AA";
			}
			else{
				mod = "NT";
			}
			var sequence:String = res[0][mod + '_sequence'];
			var type = res[0]['type'];
			if(sequence != null){
				var thisFormat:TextFormat = new TextFormat();
				var thisFont:Font = new Font1();
				thisFormat.font = thisFont.fontName;
				thisFormat.size = 12;
				thisFormat.color = 0x000000;
				
			var seqText:TextField = new TextField();
			/*verticalScroll.direction = "vertical";
			verticalScroll.setSize(seqText.width, seqText.height);
			verticalScroll.move(seqText.x + seqText.width, seqText.y);
			verticalScroll.scrollTarget = seqText;
			*/
				seqText.defaultTextFormat = thisFormat;
				seqText.embedFonts = true;
				seqText.antiAliasType = AntiAliasType.ADVANCED;
			seqText.text = ">" + res[0]['description'] + "\n\n" + sequence;
			seqText.y = space;
			seqText.width = detailPanelWidth - 60;
			seqText.multiline = true;
			seqText.wordWrap = true;
			seqText.autoSize = TextFieldAutoSize.LEFT;
			seqText.name = 'seqText';
			totalHeight = seqText.height + space;
			detailContainer.addChild(seqText);
			//update the scroll pane
			seqDetailPane.update();
			}
			else{
				if(this.seqType == "AA"){
					//get the NT sequence
					var responder = new Responder(startTranslation);
					var connection = new NetConnection();
					connection.connect(MovieClip(root).gateway);
					var thisQuery = "SELECT NT_sequence FROM " + this.projectID + "_sequences WHERE seq_id= '" + this.seqID + "'";
					connection.call("DBI.execSelect", responder, thisQuery);
				}
				else if(this.seqType == "NT"){
					
				}
			}
			}
		}
		
		private function startTranslation(res:Array){
			var responder = new Responder(showTranslatedSeqs);
			var connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var param:Array = new Array();
			param[0] = this.desc;
			param[1] = res[0]['NT_sequence'];
			connection.call("Util.translateSequence", responder, param);
		}
		
		private function showTranslatedSeqs(res:Array){
			var identifier:String = res[res.length - 1];
			var seqText:TextField = new TextField();
			var seqBackground:MovieClip = new MovieClip();
			var idText:TextField = new TextField();
			var thisFormat:TextFormat = new TextFormat();
			var thisFont:Font = new Font1();
			thisFormat.font = thisFont.fontName;
			thisFormat.size = 12;
			thisFormat.color = 0x000000;
			seqText.defaultTextFormat = thisFormat;
			seqText.embedFonts = true;
			seqText.antiAliasType = AntiAliasType.ADVANCED;
			seqText.autoSize = TextFieldAutoSize.LEFT;
			seqText.width = detailPanelWidth - 60;
			seqText.multiline = true;
			seqText.wordWrap = true;
			seqText.name = 'seqText';
			var seen:Number = 1;
			var mod:String;
			for(var i = 0; i < res.length - 1; i++){
				if(String(res[i]).indexOf(">") == 0){
					if(seen == 1){
						mod = "frame 1";
					}
					else if(seen == 2){
						mod = "frame 2";
					}
					else if(seen == 3){
						mod = "frame 3";
					}
					else if(seen == 4){
						mod = "frame -1";
					}
					else if(seen == 5){
						mod = "frame -2";
					}
					else if(seen == 6){
						mod = "frame -3";
					}
					seen++;
					seqText.appendText(">" + identifier + " " + mod + "\n");
				}
				else{
					seqText.appendText(res[i] + "\n");
				}
			}
			detailContainer.addChild(seqText);
			seqDetailPane.update();
		}
		
		private function showSeqOptions(res:Array){
		}	
		
		private function sendToBlast(e:MouseEvent):void{
			var blastPanel:BlastPanel = BlastPanel(this.parent.getChildByName('blast_panel'));
			var seqText:TextField = TextField(detailContainer.getChildByName('seqText'));
			var seq:String = seqText.text;
			var seqArray:Array = seq.split(/\r+/);
			blastPanel.addQuerySeq(seqArray[0], seqArray[1]);			
		}
		
		private function changeType(e:MouseEvent){
			if(this.seqType == 'NT'){
				this.seqType = 'AA';
			}
			else{
				this.seqType = 'NT';
			}
			this.doTypeChange(this.projectID, this.desc, this.seqType);
		}
		
		private function doTypeChange(projectID:Number, desc:String, type:String){
			if(this.getChildByName('changeText')){
				this.removeChild(this.getChildByName('changeText'));
			}
			var responder:Responder = new Responder(showSeq);
			var connection:NetConnection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query = "SELECT seq_id, description, sequence, type FROM " + projectID + "_sequences WHERE description = '" + desc + "' and type='" + type + "'";
			connection.call("DBI.execSelect", responder, query);
		}
		
		/*private function checkScrollBar():void{
			if(detailContainer.height > detailMask.height){
				scrollBar = new VScrollBar(this.detailContainer, detailContainer.height, detailMask.height);
				scrollBar.x = detailContainer.width + 1;
				scrollBar.y = detailContainer.y;
				this.addChild(scrollBar);
			}
			else{
				if(this.contains(scrollBar)){
					this.removeChild(scrollBar);
				}
			}
		}*/
		
		private function clearPanel():void{			
		}
		
	}
	
	
	
}