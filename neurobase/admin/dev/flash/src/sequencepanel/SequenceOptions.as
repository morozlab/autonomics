package sequencepanel{
	
	import flash.display.Sprite;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.events.MouseEvent;
	import flash.events.Event;
	import flash.display.DisplayObject;
	
	public class SequenceOptions extends Sprite{
		
		private var viewSeq:TextClip;
		private var addSeq:TextClip;
		private var mySequence:Sequence;
		
		public function SequenceOptions(rent:DisplayObject, sequence:Sequence){
			this.mySequence = sequence;
			this.viewSeq = new TextClip('View Sequence');
			this.viewSeq.x = 0;
			this.viewSeq.y = 0;
			viewSeq.addEventListener(MouseEvent.CLICK, sequence.displayDetails);
			this.addChild(viewSeq);
			this.addSeq = new TextClip('BLAST Sequence');
			this.addSeq.x = viewSeq.x + viewSeq.width + 5;
			this.addSeq.y = viewSeq.y;
			addSeq.addEventListener(MouseEvent.CLICK, sequence.addToBlast);
			this.addChild(addSeq);
			rent.addEventListener(ExpandEvent.EXPAND_EVENT, onExpand);
		}
		
		private function onExpand(e:ExpandEvent){
			if(mySequence.seqID > e.ID){
				this.y += e.amount;
			}
		}
		
	}
	
}