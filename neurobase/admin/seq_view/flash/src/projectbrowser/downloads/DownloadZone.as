package projectbrowser.downloads{
	
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.display.Sprite;
	import flash.events.MouseEvent;
	import flash.net.URLRequest;
	import flash.net.URLLoader;
	import flash.net.URLLoaderDataFormat;
	import flash.net.FileReference;
	import flash.events.Event;
	import flash.utils.ByteArray;
	import flash.net.NetConnection;
	import flash.net.Responder;
	
	public class DownloadZone extends Sprite{
		
		public var statText = new TextField();
		private var dlEnabled = true;
		private var binLoader:URLLoader;
		private var fileRef = new FileReference();
		private var connection;
		private var responder;
		private var compType;
		public var dlPath:String;
		public var truePath:String;
		
		public function DownloadZone():void{
			//var downloadButton:DownloadButton = new DownloadButton();
			//downloadButton.x = (300 - downloadButton.width)/2;
			//addChild(downloadButton);
			this.name = 'dlZone';
		}
		
		public function setText(str:String):void{
			statText.text = str;
		}
		
		public function enable():void{
			dlEnabled = true;
		}

		public function disable():void{
			dlEnabled = false;
		}
		
		public function checkEnabled():Boolean{
			return dlEnabled;
		}
		
		//presents the file zipped by compressFile for download
		public function presentFile(zipped:String):void{
			//load zipped file as a binary
			var splitZipped = zipped.split("/");
			truePath = swapEndPath(splitZipped[splitZipped.length - 1], dlPath);
			binLoader = new URLLoader();
			binLoader.dataFormat = URLLoaderDataFormat.BINARY;
			var binFile = new URLRequest(zipped);
			binLoader.load(binFile);
			binLoader.addEventListener(Event.COMPLETE, addStartDL);
			
		}
		
		private function swapEndPath(newPathEnd:String, oldPath:String):String{
			var splitOld = oldPath.split("\\");
			var retPath:String = "";
			for(var i = 0; i < splitOld.length; i++){
				if(splitOld[i - 1] == 'seq_view'){
					retPath = retPath +  "tmp\\" + newPathEnd;
					return retPath;
				}
				else{
					retPath = retPath + splitOld[i] + "\\";
				}
			}
			return null;
		}
		
		private function addStartDL(e:Event):void{
			var splitPath = dlPath.split("\\");
			var fileName = splitPath[splitPath.length - 1] + ".zip";
			this.setText(fileName + " ready for download! Click to start.");
			//this.addEventListener(MouseEvent.MOUSE_DOWN, startDL);
			this.startDL(new MouseEvent(MouseEvent.CLICK));
		}
		
		private function startDL(e:MouseEvent):void{
			this.removeEventListener(MouseEvent.MOUSE_DOWN, startDL);
			this.enable();
			var saveFile:ByteArray = binLoader.data;
			var splitPath = dlPath.split("\\");
			var fileName = splitPath[splitPath.length - 1] + ".zip";
			this.setText(fileName + " currently downloading!");
			fileRef.addEventListener(Event.COMPLETE, handleCompletion);
			fileRef.addEventListener(Event.CANCEL, handleCompletion);
			fileRef.save(saveFile, fileName);
		}
	
		private function handleCompletion(e:Event):void{
			//remove event listeners
			if(e.type == "cancel"){
				compType = "Cancel";
			}
			else{
				compType = "Complete";
			}
			e.currentTarget.removeEventListener(Event.COMPLETE, handleCompletion);
			e.currentTarget.removeEventListener(Event.CANCEL, handleCompletion);
			//delete temporary file from server
			responder = new Responder(deleteSuccess);
			connection = new NetConnection();
			connection.connect('http://localhost/amfphp/gateway.php');
			connection.call('FlashDownload.removeTemp', responder, truePath);
		}
		
		private function deleteSuccess(str:String):void{
			if(str == "success"){
				//enable another download
				if(compType == "Complete"){
					this.setText("Download complete. Drag an item here to download it!");
				}
				else{
					this.setText("Drag an item here to download it!");
				}
				this.enable();
			}
			else{
				//tell them to contact me
				this.setText("Error, unable to remove temporary project file. Please contact site administrator.");
			}
		}
		
	}
	
}