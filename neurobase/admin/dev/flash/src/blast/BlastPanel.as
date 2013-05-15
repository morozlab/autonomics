package blast{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFieldType;
	import flash.text.TextFormat;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	import sequencepanel.*;
	import util.*;
	import flash.events.MouseEvent;
	import flash.events.Event;
	import flash.events.DataEvent;
	import fl.controls.ComboBox;
	import fl.controls.CheckBox;
	import flash.net.Responder;
	import flash.net.NetConnection;
	import flash.net.FileReference;
	import flash.net.FileFilter;
	import flash.net.URLRequest;
	import flash.net.URLRequestMethod;
	import fl.transitions.Tween;
	import fl.transitions.easing.*; 
	import fl.transitions.TweenEvent;
	import flash.utils.Timer;
	import detail.*;
		
	public class BlastPanel extends MovieClip{
		
		private var blastTitle:BlastHeader = new BlastHeader();
		private var pasteSeqs:TextField = new TextField();
		private var blastWidth = 253;
		private var blastHeight = 300;
		private var leftPadding = 15;
		private var currentCount = 0;
		private var blastLimit = 100;
		private var connection:NetConnection;
		private var responder:Responder;
		private var projList:ComboBox = new ComboBox();
		private var blastList:ComboBox = new ComboBox();
		private var eInput:ComboBox = new ComboBox();
		private var fileUploaded:Boolean = false;
		private var serverFileName:String;
		private var manualSeqCount = 0;
		private var cover:MovieClip = new MovieClip();
		private var coverText:TextField = new TextField();
		private var initX:Number = -1;
		public var seqIndex:Number = 0;
		private var fr:FileReference = new FileReference();
		private var seqs:Array = new Array();
		private var viewText = new TextField();
		private var qTimer:Timer = new Timer(1000);
		private var qText:TextField = new TextField();
		private var statText:TextField = new TextField();
		private var statButton:MovieClip = new MovieClip();
		private var goButton:GoButton= new GoButton();
		private var bFooter:BlastFooter = new BlastFooter();
		private var curBlastID:Number = 0;
		private var blastDisplay:BlastDisplay = new BlastDisplay(800);
		private var pasteText:TextField = new TextField();
		private var useCFasta:Boolean = false;
		private var jobID:Number;
		private var uploadLimit:Number = 1024 * 1024 * 10; //10MB
		public var fileCheck:CheckBox = new CheckBox();
		public var customCheck:CheckBox = new CheckBox();
	
		public function BlastPanel():void{
		
		}
		
		public function finishUI(){
			var blastBg:MovieClip = new MovieClip();
			blastBg.graphics.beginFill(0x000000, 0);
			/*blastBg.graphics.drawRect(0, 0, blastWidth, 1);
			blastBg.graphics.drawRect(0, 0, 1, blastHeight);
			blastBg.graphics.drawRect(0, blastHeight, blastWidth, 1);
			blastBg.graphics.drawRect(blastWidth, 0, 1, blastHeight);
			*/
			blastBg.graphics.drawRect(0, 0, blastWidth, blastHeight);
			blastBg.graphics.endFill();
			this.addChild(blastBg);
			this.name = "blast_panel";
			this.addChild(blastTitle);
			var tmp:MovieClip = MovieClip(this.getChildByName('dbText'));
			var blastLabel:TextField = new TextField();
			var firstFont:Font = new Font1();
			var firstFormat:TextFormat = new TextFormat();
			firstFormat.font = firstFont.fontName;
			firstFormat.size = 14;
			firstFormat.color = 0x385F77;
			blastLabel.defaultTextFormat = firstFormat;
			blastLabel.embedFonts = true;
			blastLabel.antiAliasType = AntiAliasType.ADVANCED;
			blastLabel.x = leftPadding;
			blastLabel.y = blastTitle.y + blastTitle.height + 10;
			blastLabel.autoSize = TextFieldAutoSize.LEFT;
			blastLabel.width = 200;
			blastLabel.name = 'blastLabel';
			blastLabel.text = "BLAST Program";
			this.addChild(blastLabel);
			//finish up the blast UI
			/*//query the database for a list of projects
			responder = new Responder(finishUI);
			connection = new NetConnection();
			connection.connect(gateway);
			var query = "SELECT project_name, projectID FROM project_directory WHERE projectID IN (SELECT DISTINCT projectID FROM project_files) ORDER BY project_name";
			connection.call("DBI.execSelect", responder, query);
			*/
			var orFormat:TextFormat = new TextFormat();
			orFormat.font = firstFont.fontName;
			orFormat.size = 14;
			//orFormat.color = 0xEA9A00;
			orFormat.color = 0x385F77;
			orFormat.bold = true;
			blastList.addItem({label:'blastn', data:'blastn'});
			blastList.addItem({label:'blastp', data:'blastp'});
			blastList.addItem({label:'blastx', data:'blastx'});
			blastList.addItem({label:'tblastn', data:'tblastn'});
			blastList.addEventListener(Event.CHANGE, reloadBlastProjects);
			blastList.x = leftPadding * 2;
			blastList.width = 100;
			blastList.y = blastLabel.y + blastLabel.height + 5;
			this.addChild(blastList);
			var dbText:TextField = new TextField();
			dbText.defaultTextFormat = firstFormat;
			dbText.embedFonts = true;
			dbText.antiAliasType = AntiAliasType.ADVANCED;
			dbText.autoSize = TextFieldAutoSize.LEFT;
			dbText.text = "BLAST Project"
			dbText.x = leftPadding;
			dbText.y = blastList.y + blastList.height + 15;
			dbText.name = 'dbText';
			this.addChild(dbText);
			//populate the project list for the default BLAST program
			populateBlastList('blastn');
			projList.x = leftPadding * 2;
			projList.y = dbText.y + dbText.height + 5;
			projList.width = 200;
			this.addChild(projList);
			var evalText:TextField = new TextField();
			evalText.defaultTextFormat = firstFormat;
			evalText.embedFonts = true;
			evalText.antiAliasType = AntiAliasType.ADVANCED;
			evalText.autoSize = TextFieldAutoSize.LEFT;
			evalText.text = "Select E-Value";
			evalText.x = leftPadding;
			evalText.y = projList.y + projList.height + 15;
			this.addChild(evalText);
			eInput.width = 100;
			eInput.x = leftPadding * 2;
			eInput.y = evalText.y + evalText.height + 5;
			eInput.name = 'evalInput';
			eInput.addItem({label:"1000", data:"1000"});
			eInput.addItem({label:"100", data:"100"});
			eInput.addItem({label:"10", data:"10"});
			eInput.addItem({label:"1", data:"1"});
			eInput.addItem({label:".01", data:".01"});
			eInput.addItem({label:".001", data:".001"});
			eInput.addItem({label:"1e-04", data:'1e-04'});
			eInput.addItem({label:"1e-10", data:'1e-10'});
			eInput.addItem({label:"1e-20", data:'1e-20'});
			eInput.addItem({label:"1e-40", data:'1e-40'});
			eInput.addItem({label:"1e-50", data:'1e-50'});
			eInput.selectedIndex = 6;
			this.addChild(eInput);
			
			/*var inputText:TextField = new TextField();
			inputText.text = "Add Query Sequences (Limit 100 sequences):";
			inputText.x = leftPadding + 5;
			inputText.y = blastList.y + blastList.height + 25;
			inputText.autoSize = TextFieldAutoSize.LEFT;
			this.addChild(inputText);
			*/
			var pasteButton:MovieClip = new MovieClip();
			pasteButton.buttonMode = true;
			pasteButton.mouseChildren = false;
			var paste:TextField = new TextField();
			paste.defaultTextFormat = firstFormat;
			paste.embedFonts = true;
			paste.antiAliasType = AntiAliasType.ADVANCED;
			paste.text = "Paste Your Sequences";
			paste.autoSize = TextFieldAutoSize.LEFT;
			paste.name = 'pasteLabel';
			pasteButton.x = leftPadding;
			pasteButton.y = eInput.y + eInput.height + 30;
			pasteButton.addChild(paste);
			this.addChild(pasteButton);
			pasteText.border = true;
			pasteText.borderColor = 0x0099CC;
			pasteText.width = 220;
			pasteText.x = leftPadding + 5;
			pasteText.height = 100;
			pasteText.wordWrap = true;
			pasteText.multiline = true;
			pasteText.type = TextFieldType.INPUT;
			pasteText.y = pasteButton.y + pasteButton.height + 5;
			this.addChild(pasteText);
			var or1:TextField = new TextField();
			or1.defaultTextFormat = orFormat;
			or1.embedFonts = true;
			or1.antiAliasType = AntiAliasType.ADVANCED;
			or1.text = "Or";
			or1.autoSize = TextFieldAutoSize.LEFT;
			or1.x = leftPadding;
			or1.y = pasteText.y + pasteText.height + 15;
			this.addChild(or1);
			var uploadText:TextField = new TextField();
			uploadText.x = or1.x + or1.width + 3;
			uploadText.y = or1.y;
			uploadText.defaultTextFormat = firstFormat;
			uploadText.embedFonts = true;
			uploadText.antiAliasType = AntiAliasType.ADVANCED;
			uploadText.text = "Use File";
			uploadText.autoSize = TextFieldAutoSize.LEFT;
			this.addChild(uploadText);
			var uploadButton:UploadButton = new UploadButton();
			uploadButton.stop();
			uploadButton.buttonMode = true;
			uploadButton.mouseChildren = false;
			uploadButton.addEventListener(MouseEvent.MOUSE_OVER, uploadButton.uploadOver);
			uploadButton.addEventListener(MouseEvent.MOUSE_OUT, uploadButton.uploadOut);
			uploadButton.addEventListener(MouseEvent.CLICK, uploadFile);
			fileCheck.x = uploadText.x + uploadText.width + 5;
			fileCheck.y = uploadText.y;
			fileCheck.label = "";
			fileCheck.addEventListener(MouseEvent.CLICK, toggleFile);
			this.addChild(fileCheck);
			/*var upload:TextField = new TextField();
			upload.htmlText = "<u>upload file</u>";
			upload.autoSize = TextFieldAutoSize.LEFT;
			uploadButton.addChild(upload);
			*/
			uploadButton.name = 'uploadButton';
			uploadButton.x = leftPadding * 2;
			uploadButton.y = uploadText.y + uploadText.height + 5;
			this.addChild(uploadButton);
			var or2:TextField = new TextField();
			or2.defaultTextFormat = orFormat;
			or2.embedFonts = true;
			or2.antiAliasType = AntiAliasType.ADVANCED;
			or2.text = "Or";
			or2.autoSize = TextFieldAutoSize.LEFT;
			or2.x = leftPadding;
			or2.y = uploadButton.y + uploadButton.height + 15;
			this.addChild(or2);
			var customText:TextField = new TextField();
			customText.defaultTextFormat = firstFormat;
			customText.embedFonts = true;
			customText.antiAliasType = AntiAliasType.ADVANCED;
			customText.autoSize = TextFieldAutoSize.LEFT;
			customText.text = "Use Custom FASTA";
			customText.x = or2.x + or2.width + 3;
			customText.y = or2.y;
			this.addChild(customText);
			customCheck.x = customText.x + customText.width + 5;
			customCheck.y = customText.y;
			customCheck.label = "";
			customCheck.addEventListener(MouseEvent.CLICK, toggleCustom);
			this.addChild(customCheck);
			var emailText:TextField = new TextField();
			emailText.defaultTextFormat = firstFormat;
			emailText.embedFonts = true;
			emailText.antiAliasType = AntiAliasType.ADVANCED;
			emailText.text = "email address (for results)";
			emailText.x = leftPadding;
			emailText.y = customText.y + customText.height + 30;
			emailText.width = 230;
			emailText.autoSize = TextFieldAutoSize.LEFT;
			emailText.wordWrap = true;
			emailText.multiline = true;
			this.addChild(emailText);
			var emailInput:TextField = new TextField();
			emailInput.type = TextFieldType.INPUT;
			emailInput.y = emailText.y + emailText.height + 3;
			emailInput.x = leftPadding * 2;
			emailInput.name = 'emailInput';
			emailInput.border = true;
			emailInput.borderColor = 0x0099CC;
			emailInput.width = 190;
			emailInput.height = 20;
			emailInput.borderColor = 0x000000;
			this.addChild(emailInput);
			goButton.x = leftPadding;
			goButton.y = emailInput.y + emailInput.height + 30;
			goButton.addEventListener(MouseEvent.CLICK, checkBlastParams);
			goButton.buttonMode = true;
			//goButton.addEventListener(MouseEvent.CLICK, displayBlastResult);
			this.addChild(goButton);
			bFooter.y = goButton.y + goButton.height + 30;
			this.addChild(bFooter);
			//initialize the cover and cover text
			cover.graphics.beginFill(0xCCCCCC, .8);
			cover.graphics.drawRect(0, 0, this.width, this.height);
			cover.graphics.endFill();
			cover.name = 'blastCover';
			cover.x = this.x;
			cover.y = this.y;
			coverText.text = "Working, please wait...";
			coverText.autoSize = TextFieldAutoSize.LEFT;
			coverText.x = (cover.width - coverText.width)/2;
			coverText.y = (cover.height - coverText.height) / 2;
			cover.addChild(coverText);
			qText.multiline = true;
			statText.multiline = true;
			qText.wordWrap = true;
			statText.wordWrap = true;
			qText.autoSize = TextFieldAutoSize.LEFT;
			statText.autoSize = TextFieldAutoSize.LEFT;
			qText.width = 200;
			statText.width = 200;
			statButton.addChild(statText);
			
		}
		
		private function toggleCustom(e:MouseEvent):void{
			if(customCheck.selected == true){
				fileCheck.selected = false;
			}
		}
		
		private function toggleFile(e:MouseEvent):void{
			if(fileCheck.selected == true){
				customCheck.selected = false;
				//start file upload
				uploadFile(new MouseEvent(MouseEvent.CLICK));
			}
		}
		
		private function reloadBlastProjects(e:Event):void{
			var tmp:ComboBox = ComboBox(e.target);
			populateBlastList(tmp.selectedLabel);
		}
		
		private function populateBlastList(prog:String):void{
			var query:String;
			if(prog == "blastn" || prog == "tblastn"){
				query = "SELECT project_name, projectID FROM project_directory WHERE num_NT_seqs != '0'";
			}
			else{
				query = "SELECT project_name, projectID FROM project_directory WHERE num_AA_seqs != '0'";
			}
			query = query + " ORDER BY project_name";
			responder = new Responder(showProjectResults);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function showProjectResults(res:Array):void{
			projList.removeAll();
			for(var i = 0; i < res.length; i++){
				projList.addItem({label:res[i]['project_name'], data:res[i]['projectID']});
			}
		}
		
		
		private function uploadFile(e:MouseEvent){
			fr.addEventListener(Event.SELECT, fileSelectHandler);
			fr.addEventListener(Event.CANCEL, cancelHandler);
			fr.browse(new Array(new FileFilter("FASTA Files", "*.fas;*.fna;*.fsa;*.fasta;*.txt)")));
		}
		
		private function cancelHandler(e:Event):void{
			fileCheck.selected = false;
			fileUploaded = false;
			fr.removeEventListener(Event.SELECT, fileSelectHandler);
			fr.removeEventListener(Event.CANCEL, cancelHandler);
		}
		
		private function fileSelectHandler(e:Event){
			fr.removeEventListener(Event.SELECT, fileSelectHandler);
			fr.removeEventListener(Event.CANCEL, cancelHandler);
			if(fr.size > uploadLimit){
				displayError('upSize');
				fileUploaded = false;
				fileCheck.selected = false;
			}
			else{
				fileUploaded = true;
				fileCheck.selected = true;
				customCheck.selected = false;
			}
		}
		
		private function uploadCompleteHandler(e:DataEvent){
			var fileName:String = e.data
			responder = new Responder(blastPrepared, blastFail);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var param:Array = new Array(5);
			var seqs:String = pasteSeqs.text;
			var pattern:RegExp = /\r/g;
			seqs = seqs.replace(pattern, "\n");
			param[0] = seqs;
			param[1] = fileName;
			param[2] = projList.selectedItem.data;
			param[3] = eInput.selectedItem.data;
			param[4] = TextField(this.getChildByName('emailInput')).text
			param[5] = blastList.selectedItem.data;
			connection.call("BLAST.prepareBLAST", responder, param);
		}
		
		private function checkBlastParams(e:MouseEvent):void{
			if(pasteText.text == "" && fileCheck.selected == false && customCheck.selected == false){
				displayError('query');
			}
			else{
				var okToContinue:Boolean = true;
				var errCode:String = "";
				//add the BLAST cover
				this.addChild(cover);
				if(customCheck.selected == true){
					//check the sequences in the custom fasta to see if they are all one or the other
					var bPanel:BottomPanel = BottomPanel(this.parent.getChildByName('bottom_panel'));
					var tmp:Array = bPanel.getCustomSeqs();
					if(tmp.length > 0){
						var seenNT:Boolean = false;
						var seenAA:Boolean = false;
						for(var i = 0; i < tmp.length; i++){
							if(tmp[i][2] == "NT"){
								seenNT = true;
								if(seenAA == true){
									okToContinue = false;
									errCode = 'cMisMatch';
									break;
								}
							}
							else{
								seenAA = true;
								if(seenNT == true){
									okToContinue = false;
									errCode = 'cMisMatch';
									break;
								}
							}					
						}
					}
					else{
						okToContinue = false;
						errCode = 'cEmpty';
					}
						
				}
				if(okToContinue == true){
						parsePastedSequences();
				}
				else{
					this.removeChild(cover);
					displayError(errCode);
				}
				/*responder = new Responder(checkDBType);
				connection = new NetConnection();
				connection.connect(inc.gateway);
				var query:String = "SELECT num_NT_seqs, num_AA_seqs FROM project_directory WHERE projectID ='" + projList.selectedItem.data + "'";
				connection.call("DBI.execSelect", responder, query);
				*/
			}
		}
		
		private function removeErrorBox(e:MouseEvent):void{
			this.removeChild(this.getChildByName('errorBox'));
		}
		
		private function checkDBType(res:Array):void{
			var blastP:String = blastList.selectedItem.data;
			if((blastP == 'blastn' || blastP == 'tblastn') && res[0]['num_NT_seqs'] == 0){
				displayError('noNT');
			}
			else if((blastP == 'blastp' || blastP == 'blastx') && res[0]['num_AA_seqs'] == 0){
				displayError('noAA');
			}
			else{
				parsePastedSequences();
			}
		}
		
		private function displayError(err:String):void{
			if(this.getChildByName('errorBox')){
				this.removeChild(this.getChildByName("errorBox"));
			}
			var errorBox:MovieClip = new MovieClip();
			errorBox.name = 'errorBox';
			errorBox.y = goButton.y - 50;
			var errFormat:TextFormat = new TextFormat();
			var errFont:Font = new Font1();
			errFormat.font = errFont.fontName;
			errFormat.size = 10;
			errFormat.color = 0xFB0000;
			var errorText:TextField = new TextField();
			errorText.x = 10;
			errorText.width = 180;
			errorText.y = 10;
			errorText.wordWrap = true;
			errorText.multiline = true;
			errorText.defaultTextFormat = errFormat;
			errorText.embedFonts = true;
			errorText.antiAliasType = AntiAliasType.ADVANCED;
			errorText.autoSize = TextFieldAutoSize.LEFT;
			errorBox.addChild(errorText);
			var errorClip:TextClip = new TextClip('OK');
			errorClip.x = (200 - errorClip.width) /2;
			errorClip.addEventListener(MouseEvent.CLICK, removeErrorBox);
			errorBox.addChild(errorClip);
			if(err == 'cMisMatch'){
				errorText.text = "ERROR: sequence mismatch in custom FASTA";
			}
			else if(err == "cEmpty"){
				errorText.text = "ERROR: custom FASTA empty";
				
			}
			else if(err == 'upSize'){
				errorText.text = "ERROR: file exceeds maximum size (10MB)";
			}
			else if(err == 'query'){
				errorText.text = "ERROR: no set of query sequences found";
			}
			else if(err == 'pProgMisMatch'){
				
			}
			errorClip.y = errorText.y + errorText.height + 3;
			errorBox.graphics.lineStyle(1, 0x000000, 1);
			errorBox.graphics.beginFill(0xFFFFFFF, 1);
			errorBox.graphics.drawRect(0, 0, 200, errorBox.height + 10);
			errorBox.graphics.endFill();
			this.addChild(errorBox);
		}
		
		private function prepareBLAST(){
			if(customCheck.selected == true){
				//get the bottom display object
				var bDisplay:BottomPanel = BottomPanel(this.parent.getChildByName('bottom_panel'));
				seqs = bDisplay.getCustomSeqs();
			}
			responder = new Responder(blastPrepared, blastFail);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var param:Array = new Array(6);
			if(fileCheck.selected == false){
				param[0] = seqs;
				param[1] = 'none';
				param[2] = projList.selectedItem.data;
				param[3] = eInput.selectedItem.data;
				param[4] = TextField(this.getChildByName('emailInput')).text
				param[5] = blastList.selectedItem.data;
				connection.call("BLAST.prepareBLAST", responder, param);
			}
			else{
				var req:URLRequest = new URLRequest(MovieClip(root).basePath + '/admin/seq_view/src/flashUpload.php');
				req.method = URLRequestMethod.POST;
				fr.addEventListener(DataEvent.UPLOAD_COMPLETE_DATA, uploadCompleteHandler);
				fr.upload(req, 'flashFile');
			}
		}
		
		private function blastPrepared(res:Array){
			//remove the BLAST cover
			this.removeChild(this.getChildByName('blastCover'));
			//remove the result button if it's there
			if(this.getChildByName('gResults')){
				this.removeChild(this.getChildByName('gResults'));
			}
			//set the cur BLASTID for queue purposes
			curBlastID = res[1];
			//set the jobID to pass to display
			jobID = res[0];
			this.stage.addEventListener(Event.ENTER_FRAME, checkTimer);
			//add text about BLAST job
			var qFormat:TextFormat = new TextFormat();
			var qFont:Font = new Font1();
			qFormat.font = qFont.fontName;
			qFormat.size = 14;
			qFormat.color = 0x000000;
			qText.defaultTextFormat = qFormat;
			qText.embedFonts = true;
			qText.antiAliasType = AntiAliasType.ADVANCED;
			qText.setTextFormat(qFormat);
			qText.htmlText = "Job " + String(jobID);
			qFormat.color = 0xE6880F;
			if(!(this.contains(qText))){
				//add the blast divider
				var divider:BlastDivider = new BlastDivider();
				divider.y = goButton.y + goButton.height + 20;
				this.addChild(divider);
				qText.x = goButton.x;
				qText.y = divider.y + divider.height + 5;
				qText.autoSize = TextFieldAutoSize.LEFT;
				this.addChild(qText);
				//make the job text
			}
			//look up where the BLAST job is in the queue
			responder = new Responder(lowestBlastID);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			var query:String = "SELECT min(blast_id)as lowest FROM blast_queue";
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function lowestBlastID(res:Array){
			if(!(this.contains(statButton))){
				statButton.x = qText.x;
				statButton.y = qText.y + qText.height;
				this.addChild(statButton);
				
			}
			var statFormat:TextFormat = new TextFormat();
			var statFont:Font = new Font1();
			statFormat.font = statFont.fontName;
			statFormat.size = 14;
			statText.embedFonts = true;
			statText.antiAliasType = AntiAliasType.ADVANCED;
			var gResults:GetResults = new GetResults();
			if(res[0]['lowest'] != null){
				if(curBlastID < res[0]['lowest']){
					this.stage.removeEventListener(Event.ENTER_FRAME, checkTimer);
					while(statButton.hasEventListener(MouseEvent.CLICK)){
						statButton.removeEventListener(MouseEvent.CLICK, displayBlastResult);
					}
					
					statText.htmlText = "";
					gResults.buttonMode = true;
					gResults.x = statButton.x;
					gResults.y = statButton.y;
					gResults.name = 'gResults';
					gResults.addEventListener(MouseEvent.CLICK, displayBlastResult);
					this.addChild(gResults);
					qTimer.reset();
				}
				else{
					while(statButton.hasEventListener(MouseEvent.CLICK)){
						statButton.removeEventListener(MouseEvent.CLICK, displayBlastResult);
					}
					statText.autoSize = TextFieldAutoSize.LEFT;
					statText.htmlText = "Position " + (curBlastID - res[0]['lowest'] + 1) + " in the queue.";
					qTimer.start();
				}
			}
			else{
				while(statButton.hasEventListener(MouseEvent.CLICK)){
						statButton.removeEventListener(MouseEvent.CLICK, displayBlastResult);
				}
				qTimer.reset();
				this.stage.removeEventListener(Event.ENTER_FRAME, checkTimer);
				statText.htmlText = "";
				gResults.buttonMode = true;
				gResults.x = statButton.x;
				gResults.y = statButton.y + 3;
				gResults.name = 'gResults';
				gResults.addEventListener(MouseEvent.CLICK, displayBlastResult);
				this.addChild(gResults);
			}
			statText.setTextFormat(statFormat);
			bFooter.y = statButton.y + 60;
		}
		
		private function displayBlastResult(e:MouseEvent):void{
			this.parent.addChild(blastDisplay);
			blastDisplay.initDisplay(jobID);
		}
		
		private function checkTimer(e:Event):void{
			if(qTimer.currentCount >= 2){
				qTimer.reset();
				//check to see if the job is still in the queue
				responder = new Responder(lowestBlastID);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				var query = "SELECT min(blast_id) as lowest FROM blast_queue";
				connection.call("DBI.execSelect", responder, query);
			}
		}
		
		private function blastFail(res:Object){
			
		}
		
		public function showPasteBox(e:MouseEvent){
			var tmp:TextField = new TextField();
			tmp.x = 10;
			tmp.y = 5;
			tmp.autoSize = TextFieldAutoSize.LEFT;
			tmp.text = "Please paste your sequences below:";
			var pasteBg:MovieClip = new MovieClip();
			pasteBg.graphics.beginFill(0xFFFFFF, 1);
			pasteBg.graphics.drawRect(0,0, 300, 600);
			pasteBg.graphics.endFill();
			/*pasteBg.graphics.beginFill(0x000000, 1);
			pasteBg.graphics.drawRect(0, 0, 400, 1);
			pasteBg.graphics.drawRect(0, 0, 1, 600);
			pasteBg.graphics.drawRect(0, 600, 400, 1);
			pasteBg.graphics.drawRect(400, 0, 1, 600);
			*/pasteSeqs.type = TextFieldType.INPUT;
			pasteSeqs.border = true;
			pasteSeqs.borderColor = 0x000000;
			pasteSeqs.width = 260;
			pasteSeqs.height = 500;
			pasteSeqs.multiline = true;
			pasteSeqs.wordWrap = true;
			pasteSeqs.x = 20;
			pasteSeqs.y = tmp.y + tmp.height + 5;
			pasteBg.addChild(tmp);
			pasteBg.addChild(pasteSeqs);
			pasteBg.x = 0
			var finishText:TextField = new TextField();
			finishText.text = 'Finished';
			finishText.autoSize = TextFieldAutoSize.LEFT;
			var finishButton = new MovieClip();
			finishButton.graphics.beginFill(0xFFFFFF, 1);
			finishButton.graphics.drawRect(0, 0, finishText.width + 10, 15);
			finishButton.graphics.endFill();
			finishButton.mouseChildren = false;
			finishButton.buttonMode = true;
			finishButton.x = pasteSeqs.x;
			finishButton.y = pasteSeqs.y + pasteSeqs.height + 10;
			finishButton.addChild(finishText);
			finishButton.addEventListener(MouseEvent.CLICK, parsePastedSequences);
			finishText.x = (finishButton.width - finishText.width) / 2
			pasteBg.addChild(finishButton);
			pasteBg.name = 'pasteBg';
			this.addChild(pasteBg);
		}
		
		public function parsePastedSequences(){
			seqs = new Array();
			//parse the sequences that were in the box
			var tmp:String = pasteText.text;
			var splitStr:Array = tmp.split(/\r+/);
			var seq:String = "";
			var seqID:String = "";
			var seenSeqs = 0;
			seqIndex = 0;
			for(var i = 0; i < splitStr.length; i++){
				var line:String = splitStr[i];
				if(seqIndex == 100){
					break;
				}
				if(line.match(/^>/)){
					//if this is the first sequence, set the ID equal to the line and clear the sequence string
					if(seenSeqs == 0){
						seqID = line;
						seq = "";
						seenSeqs = 1;
					}
					//this is an additional sequence, set the old info into
					else{
						seqs[seqIndex] = [seqID, seq];
						seq = "";
						seqID = line;
						seqIndex++;
					}
				}
				else{
					seq = seq + line;
				}
			}
			if(seqID.length > 0){
				seqs[seqIndex] = [seqID, seq];			
				seqIndex++;
			}
			else if(seq.length > 0 && seenSeqs == 0){
				seqs[seqIndex] = ['userseq1', seq];
				seqIndex++;
			}
			pasteSeqs.text = "";
			this.prepareBLAST();
		}
		
		public function useCustomFasta(toUse:Array):void{
			this.useCFasta = true;
		}
		
		public function addQuerySeq(ID:String, seq:String){
			if(seqIndex < blastLimit){
				//check if this sequence already exists in the array
				for(var i:Number = 0; i < seqs.length; i++){
					if(seqs[i][0] == ID){
						return;
					}
				}
				seqs[seqIndex] = [ID, seq];
				seqIndex++;
				updateViewText();
			}
		}
		
		private function updateViewText(){
			viewText.htmlText =  "<u>view (" + seqIndex + "/100)</u>";
		}
		
		private function displayQuerySeqs(e:MouseEvent){
			if(seqIndex > 0){
				//roll out the BLAST panel
				if(initX == -1){
					initX = this.x;
				}
				var tmp:Tween = new Tween(this, 'x', Regular.easeOut, this.x, this.stage.stageWidth + 10, .5, true);
				//display the sequence IDs that have been manually submitted as query sequences
				tmp.addEventListener(TweenEvent.MOTION_FINISH, doQueryDisplay);
			}
		}
		
		private function doQueryDisplay(e:TweenEvent){
			//add a new container to hold all of the sequences
			var tmpContainer:MovieClip = new MovieClip();
			tmpContainer.x = this.initX;
			tmpContainer.name = 'tmpContainer';
			var lastY:Number = 30;
			for(var i = 0; i < this.seqIndex; i++ ){
				var paste:PastedSequence = new PastedSequence(i, seqs[i][0], this.width - 40, this);
				paste.x = 0;
				paste.y = lastY + 3;
				lastY = lastY + 3 + paste.height;
				tmpContainer.addChild(paste);
			}
			var returnButton:TextClip = new TextClip("<- Return to BLAST Panel");
			returnButton.x = this.initX;
			returnButton.y = 0;
			returnButton.name = 'returnButton';
			tmpContainer.y = returnButton.y + returnButton.height + 3;
			this.parent.addChild(returnButton);
			returnButton.addEventListener(MouseEvent.CLICK, returnToBlastPanel);
			this.parent.addChild(tmpContainer);
		}
		
		public function removePastedSequence(e:MouseEvent){
			var pasted:PastedSequence = PastedSequence(TextClip(e.currentTarget).datum);
			var tmpContainer:MovieClip = MovieClip(this.parent.getChildByName('tmpContainer'));
			//shift the removed sequence
			for(var i = pasted.thisIndex; i < seqs.length - 1; i++){
				seqs[i] = seqs[i + 1];
			}
			seqs.pop();
			//decrement the index;
			seqIndex--;
			//visually remove the sequence
			tmpContainer.removeChild(pasted);
			//update the text for how many sequences there are
			updateViewText();
		}
		
		private function returnToBlastPanel(e:MouseEvent){
			this.parent.removeChild(this.parent.getChildByName('tmpContainer'));
			this.parent.removeChild(this.parent.getChildByName('returnButton'));
			//start the tween of the blastpanel back on screen
			var tmp:Tween = new Tween(this, 'x', Regular.easeOut, this.x, this.initX, 1, true);
		}
		
	}
	
}