package sequencepanel {
	import flash.events.Event;

	public class ExpandEvent extends Event {
		
		public static  const EXPAND_EVENT:String = "expandevent";

		public function ExpandEvent(type:String, ID:int, amount:int ) {
			super(type, true);
			
			this.amount = amount
			this.ID = ID;
		}
		
		public var amount:int;
		public var ID:int;
		
		override public function clone():Event {
			return new ExpandEvent(this.type, this.ID, this.amount);
		}
	}
}