package projectbrowser{
	
	import flash.display.Sprite;
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFormat;
	import flash.events.MouseEvent;
	import fl.controls.ComboBox;
	import sequencepanel.ExpandEvent;
	
	public class CFastaSeq extends Sprite{
		
		private var trashCan:TrashCan = new TrashCan();
		private var seqContainer:Sprite = new Sprite();
		private var seqFormat:TextFormat = new TextFormat();
		private var seqFont:TahomaFont = new TahomaFont;
		private var seqText:TextField = new TextField();
		public var seqTextContainer:MovieClip = new MovieClip();
		public var seqID:Number;
		public var projID:Number;
		public var desc:String;
		private var seqType:String;
		public var addID:Number;
		private var rent:MovieClip;
		private var customFasta:CustomFasta;
		
		public function CFastaSeq(projID:Number, seqID:Number, desc:String, rent:MovieClip, addID:Number, cFasta:CustomFasta):void{
			this.projID = projID;
			this.seqID = seqID;
			this.desc = desc;
			this.rent = rent;
			this.customFasta = cFasta;
			this.addID = addID;
			this.name = String(addID);
			seqFormat.font = seqFont.fontName;
			seqFormat.size = 12;
			trashCan.buttonMode = true;
			this.addChild(trashCan);
			trashCan.addEventListener(MouseEvent.CLICK, removeThis);
			this.customFasta.addEventListener(ExpandEvent.EXPAND_EVENT, moveUp);
			seqTextContainer.x = trashCan.width + 60;
			seqTextContainer.addChild(seqText);
			seqTextContainer.mouseChildren = false;
			seqTextContainer.buttonMode = true;
			seqText.defaultTextFormat = seqFormat;
			seqText.embedFonts = false;
			seqText.y = 0;
			seqText.width = 400;
			seqText.autoSize = TextFieldAutoSize.LEFT;
			seqText.wordWrap = true;
			seqText.multiline = true;
			seqText.text = desc;
			this.addChild(seqTextContainer);
			if(seqText.height - trashCan.height > 0){
				trashCan.y = 0;
			}
			else{
				trashCan.y = (seqText.height - trashCan.height)/2;
			}
		}
		
		private function moveUp(e:ExpandEvent):void{
			if(this.addID > e.ID){
				var box:CustomComboBox = CustomComboBox(rent.getChildByName("comboBox" + this.addID));
				this.addID -= 1;
				//change the ID for the combo box
				box.customSeqID = this.addID;
				box.name = "comboBox" + this.addID;
				this.name = String(addID);
				this.y -= e.amount;
				//move the combo box for this sequence
				box.y = this.y;
			}
		}
		
		private function removeThis(e:MouseEvent):void{
			this.customFasta.removeEventListener(ExpandEvent.EXPAND_EVENT, moveUp);
			customFasta.removeSequence(this.addID);
		}
		
	}
	
	
}