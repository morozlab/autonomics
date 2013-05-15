package detail{
	
	import flash.display.Sprite;
	import flash.text.TextFormat;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.Font;
	import flash.events.MouseEvent;
	import projectbrowser.*;
	
	public class BottomPanel extends Sprite{
		
		private var detailPanel:DetailPanel = new DetailPanel();
		private var cFasta:CustomFasta = new CustomFasta();
		private var header:DetailHeader = new DetailHeader();
		private var footer:SequenceFooter = new SequenceFooter();
		private var headerText:TextField = new TextField();
		private var thisFormat:TextFormat = new TextFormat;
		private var detailLink:DetailLink = new DetailLink();
		
		public function BottomPanel():void{
			detailPanel.visible = false;
			cFasta.visible = false;
			header.gotoAndStop(1);
			//add the header
			this.addChild(header);
			detailPanel.y = header.height + 40;
			cFasta.y = header.height + 10;
			//add the two elements
			this.addChild(detailPanel);
			this.addChild(cFasta);
			//add the navigation clips
			detailLink.x = header.width - detailLink.width - 5;
			detailLink.y = cFasta.y;
			detailLink.setText("Custom FASTA");
			detailLink.mouseChildren = false;
			detailLink.buttonMode = true;
			detailLink.addEventListener(MouseEvent.MOUSE_OVER, detailLink.onMouseOver);
			detailLink.addEventListener(MouseEvent.ROLL_OUT, detailLink.onMouseOut);
			detailLink.addEventListener(MouseEvent.CLICK, toggleDisplay);
			this.addChild(detailLink);
			//add the footer
			footer.y = detailPanel.y + detailPanel.height + 10;
			this.addChild(footer);
			//draw the mask for the bottom panel
			var bottomMask:Sprite = new Sprite();
			bottomMask.graphics.beginFill(0, 0);
			bottomMask.graphics.drawRect(0, 0, this.width + 10, this.height);
			bottomMask.graphics.endFill();
			this.addChild(bottomMask);
			this.mask = bottomMask;
			this.name = 'bottom_panel';
		}
		
		private function toggleDisplay(e:MouseEvent):void{
			if(cFasta.visible == false){
				detailPanel.visible = false;
				cFasta.visible = true;
				header.gotoAndStop(2);
				detailLink.setText("<u>Seq Detail</u>");
			}
			else{
				cFasta.visible = false;
				detailPanel.visible = true;
				header.gotoAndStop(1);
				detailLink.setText("<u>Custom FASTA</u>");
			}
		}
		
		public function displaySequenceDetail(seqID:Number, projectID:Number, description:String, type:String):void{
			cFasta.visible = false;
			detailPanel.visible = true;
			header.gotoAndStop(1);
			detailLink.setText("Custom FASTA");
			this.detailPanel.displaySequenceDetail(seqID, projectID, description, type);
		}
	
		public function addCustomFastaSequence(seqID:Number, projectID:Number, desc:String, defaultType:String):void{
			detailPanel.visible = false;
			cFasta.visible = true;
			header.gotoAndStop(2);
			detailLink.setText("Seq Detail");
			this.cFasta.addSequence(seqID, projectID, desc, defaultType);
		}
		
		public function getCustomSeqs():Array{
			return cFasta.getSeqs();
		}
		
	}
	
	
	
}