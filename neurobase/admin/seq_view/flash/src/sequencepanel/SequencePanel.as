package sequencepanel{
	
	import flash.display.MovieClip;
	import flash.text.TextField;
	import flash.text.TextFieldAutoSize;
	import flash.text.TextFieldType;
	import flash.text.AntiAliasType;
	import flash.text.Font;
	import flash.events.*;
	import flash.text.TextFormat;
	import flash.text.Font;
	import flash.net.NetConnection;
	import flash.net.Responder;
	import flash.utils.*;
	import flash.geom.Rectangle;
	import flash.display.Sprite;
	import fl.containers.ScrollPaneNoBorder;
	import util.*;
	
	
	public class SequencePanel extends MovieClip{
		
		private var connection:NetConnection;
		private var responder:Responder;
		private var query:String;
		private var pageSize = 100;
		private var curPage = 1;
		private var curID;
		private var curType;
		private var lastY = 0;
		internal var lastID;
		internal var expandAmount = 30;
		private var lastPage;
		private var searchString:String;
		private var projectQuery:String;
		private var baseQuery;
		private var panelWidth = 530;
		private var panelHeight = 400;
		private var beforeSearch:String = "";
		private var maskClip:MovieClip = new MovieClip();
		private var scrollBar:MovieClip;
		private var container:MovieClip = new MovieClip();
		private var scrollRatio:Number;
		private var curScrollPosition:Number = 0;
		private var sequenceBox:MovieClip;
		private var spacer:Number = 5;
		private var loadType:String;
		private var curFileID:Number = -1;
		private var sequenceFooter:SequenceFooter = new SequenceFooter();
		private var seqPane:ScrollPaneNoBorder;
		private var searchLink:SearchLink;
		private var fromSearch:Boolean = false;
		private var thumbPath:String = "none";
		
		public function SequencePanel(){
			this.name = 'seqPanel';
		}
		
		public function initSequencePanel(){
			var myFont = new eVerdana();
			var myFormat:TextFormat = new TextFormat();
			myFormat.font = myFont.fontName;
			myFormat.size = 10;
			var hitBox = new MovieClip();
			hitBox.x = 0;
			hitBox.y = 0;
			hitBox.name = 'hitBox';
			container.addChild(hitBox);
			hitBox.graphics.beginFill(0, 0);
			hitBox.graphics.drawRect(0, 0, panelWidth, panelHeight);
			hitBox.graphics.endFill();
			hitBox.width = panelWidth;
			var seqHeader:SequenceHeader = new SequenceHeader();
			this.addChild(seqHeader);
			var pageBox = new MovieClip();
			pageBox.name = 'pageBox';
			pageBox.x = 0;
			container.addChild(pageBox);
			sequenceBox = new MovieClip();
			sequenceBox.name = 'sequenceBox';
			seqPane = new ScrollPaneNoBorder();
			seqPane.source = sequenceBox;
			seqPane.setSize(panelWidth + 20, panelHeight - 30);
			seqPane.update();
			seqPane.y = 40;
			seqPane.x = 0;
			container.addChild(seqPane);
			container.y = seqHeader.height + 35;
			sequenceFooter.y = container.y + container.height + 20;
			this.addChild(sequenceFooter);
			this.addChild(container);
		}
		
		public function setThumbPath(path:String):void{
			this.thumbPath = path;
		}
		
		private function checkForScroll(){
			if(container.height > panelHeight){
				var tmp:MovieClip = MovieClip(scrollBar.getChildByName('scrollPane'));
				tmp.y = 0
				Math.round(scrollRatio = (lastY - maskClip.height)/(scrollBar.height - 40));
				scrollBar.visible = true;
				curScrollPosition = 0;
			}
			else{
				scrollBar.visible = false;
			}
		}
		
		private function startScroll(e:MouseEvent){
			MovieClip(scrollBar.getChildByName('scrollPane')).startDrag(false, new Rectangle(0, 0, 0, scrollBar.height - 40));
			this.stage.addEventListener(MouseEvent.MOUSE_UP, endScroll);
			this.stage.addEventListener(MouseEvent.MOUSE_MOVE, doScroll);
		}
		
		private function endScroll(e:MouseEvent){
			MovieClip(scrollBar.getChildByName('scrollPane')).stopDrag();
			this.stage.removeEventListener(MouseEvent.MOUSE_UP, endScroll);
			this.stage.removeEventListener(MouseEvent.MOUSE_MOVE, doScroll);
		}
		
		internal function moveBar(){
			
		}
		
		internal function resizeScroll(opt:String, exp:Number){
			var tmp = MovieClip(scrollBar.getChildByName('scrollPane'));
			//determine the new scroll bar ratio
			scrollRatio = (sequenceBox.height - maskClip.height)/(scrollBar.height - 40);
			//figure out where the sequence box is in relation to the mask
			var dy:Number;
			if(opt == 'expand'){
				dy = Math.abs(sequenceBox.y - maskClip.y);
			   if(tmp.y != 0){
				   tmp.y = dy / scrollRatio;
			   }
			}
			else{
				if(sequenceBox.y + exp > 30){
					sequenceBox.y = 30;
				}
				else{
					sequenceBox.y += exp;
				}
				dy = Math.abs(sequenceBox.y - maskClip.y);
				if((dy / scrollRatio) > (scrollBar.height - tmp.height)){
					tmp.y = scrollBar.height - tmp.height;
				}
				else{
					tmp.y = dy / scrollRatio;
				}
			}
			//trace(maskClip.y + "|" + sequenceBox.y);
			curScrollPosition = tmp.y;
		}
		
		private function doScroll(e:MouseEvent){
			var tmp =  MovieClip(scrollBar.getChildByName('scrollPane'));
			var posChange = curScrollPosition -  tmp.y;
			curScrollPosition = tmp.y;
			sequenceBox.y = sequenceBox.y + (posChange * scrollRatio);
		}
		
		public function loadSeqs(ID:int, type:String, page:Number){
			curPage = page;
			curID = ID;
			curType = type;
			clearChildren('sequenceBox');
			clearChildren('pageBox');
			fromSearch = false;
			if(this.getChildByName('origClip')){
				this.removeChild(this.getChildByName('origClip'));
			}
			if(type == 'project'){
				responder = new Responder(handleProjectQuery);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				query = "SELECT num_AA_seqs, num_NT_seqs, projectID, project_name, default_type FROM project_directory WHERE projectID ='" + ID + "'";
				connection.call("DBI.execSelect", responder, query);				
				baseQuery = query;
			}
			else if(type == 'file'){
				//get the project name of the file in question
				responder = new Responder(handleFileQuery);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				query = "SELECT num_seqs, fileID, projectID, file_name FROM project_files WHERE fileID ='" + ID + "'";
				connection.call("DBI.execSelect", responder, query);
				
			}
			
		}
		
		private function handleProjectQuery(res:Array){
			if(res.length > 0){
				if(res[0]['num_AA_seqs'] > 0 || res[0]['num_NT_seqs'] > 0){
					//show descriptor of the project being shown
					addDescriptionText(res[0]['project_name'], 'project');
					var numAA:Number = res[0]['num_AA_seqs'];
					var numNT:Number = res[0]['num_NT_seqs'];
					var numDisplaySeqs:Number;
					if(res[0]['default_type'] == 'AA'){
						numDisplaySeqs = numAA;
						loadType = "AA";
					}
					else{
						numDisplaySeqs = numNT;
						loadType = "NT";
					}
					curPage = 1;
					lastPage = Math.ceil(numDisplaySeqs/pageSize);
					baseQuery = "SELECT seq_id, project_id, description, type FROM " + res[0]['projectID'] + "_sequences";
					var limiter:String = "limit " + ((curPage - 1) * pageSize) + ", " + pageSize; 
					query = baseQuery + " " + limiter;
					responder = new Responder(doLoad);
					connection = new NetConnection();
					connection.connect(MovieClip(root).gateway);
					connection.call("DBI.execSelect", responder, query);
				}
			}
		}
		
		private function handleFileQuery(res:Array){
			if(res.length > 0){
				//display the identifier for the file
				addDescriptionText(res[0]['file_name'], 'file');
				curID = res[0]['projectID'];
				curFileID = res[0]['fileID'];
				curPage = 1;
				lastPage = Math.ceil(res[0]['num_seqs']/pageSize);
				baseQuery = "SELECT seq_id, project_id, description, type FROM " + res[0]['projectID'] + "_sequences WHERE file_id ='" + res[0]['fileID'] + "'";
				var limiter:String = "limit " + ((curPage - 1) * pageSize) + ", " + pageSize; 
				query = baseQuery + " " + limiter;
				fromSearch = false;
				responder = new Responder(doLoad);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				connection.call("DBI.execSelect", responder, query);
			}
		}
		
		public function loadPage(ID:int, type:String, page:Number){
			curPage = page;
			clearChildren('sequenceBox');
			if(curPage < 1){
				curPage = 1;
			}
			else if(curPage > lastPage){
				curPage = lastPage;
			}
			TextClip(MovieClip(container.getChildByName('pageBox')).getChildByName('pageClip')).thisText.text = 'Page ' + curPage + ' of ' + lastPage;
			//reposition page navigation elements
			reposPageNav();	
			var limiter:String = "limit " + ((curPage - 1) * pageSize) + ", " + pageSize; 
			query = baseQuery + " " + limiter;
			responder = new Responder(doLoad);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			connection.call("DBI.execSelect", responder, query);
		}
		
		private function reposPageNav():void{
			var pageClip:TextClip = TextClip(MovieClip(container.getChildByName('pageBox')).getChildByName('pageClip'));
			var nextClip:TextClip = TextClip(MovieClip(container.getChildByName('pageBox')).getChildByName('nextClip'));
			var prevClip:TextClip = TextClip(MovieClip(container.getChildByName('pageBox')).getChildByName('prevClip'));
			//var pageInput:TextField = TextField(MovieClip(container.getChildByName('pageBox')).getChildByName('pageInput'));
			prevClip.x = pageClip.x + pageClip.width + spacer;
			//pageInput.x = prevClip.x + prevClip.width + spacer;
			nextClip.x = prevClip.x + prevClip.width + spacer;
		}
				
		private function handleCount(res:Array){
			var numRows = res.length;
			lastPage = Math.ceil(numRows/pageSize);
			if(lastPage == 0){
				lastPage = 1;
			}
			if(curPage < 1){
				curPage = 1;
			}
			else if(curPage > lastPage){
				curPage = lastPage;
			}
			var limiter:String = "limit " + ((curPage - 1) * pageSize) + ", " + pageSize; 
			query = baseQuery + " " + limiter;
			responder = new Responder(doLoad);
			connection = new NetConnection();
			connection.connect(MovieClip(root).gateway);
			connection.call("DBI.execSelect", responder, query);
			
		}
		
		private function doLoad(res:Array){
			var seqNumber = ((curPage - 1) * pageSize) + 1; 
			//remove the loading text, if it's there
			if(container.getChildByName('ldMovie')){
				container.removeChild(container.getChildByName('ldMovie'));
			}
			if(!(MovieClip(container.getChildByName('pageBox')).getChildByName('pageClip'))){
					addNavigation();
			}
			if(res.length > 0){
				//add navigation elements
				loadType = res[0]['type'];
			}
			for(var i = 0; i < res.length; i++){
				addSequence(seqNumber, res[i]['seq_id'], res[i]['project_id'], res[i]['description'], res[i]['type']);
				seqNumber++;
			}
			var fillerMovie:MovieClip = new MovieClip();
			fillerMovie.graphics.beginFill(0, 0);
			fillerMovie.graphics.drawRect(0, 0, 10, 30);
			fillerMovie.graphics.endFill();
			fillerMovie.y = lastY + 10;
			sequenceBox.addChild(fillerMovie);
			seqPane.update();
		}
		
		private function addSequence(seqNumber:Number, ID:Number, project:Number, description:String, type:String){
			var seq:Sequence = new Sequence(seqNumber, ID, project, description, this, type);
			seq.y = lastY + 20;
			seq.initYs();
			lastY = lastY + 20 + seq.height;
			sequenceBox.addChild(seq);
			lastID = seq.seqID;
		}
		
		private function addNavigation(){
			var pageBox = MovieClip(container.getChildByName('pageBox'));
			var curWidth = 0;
			//space between navigation elements
			//add the page number display
			var pageClip = new TextClip('Page ' + curPage + ' of ' + lastPage);
			pageClip.name = 'pageClip';
			pageClip.buttonMode = false;
			pageBox.addChild(pageClip);
			curWidth += pageClip.width;
			//add the previous button display
			var prevClip = new TextClip('< previous');
			prevClip.x = curWidth + spacer;
			prevClip.name = 'prevClip';
			curWidth += prevClip.width + spacer;
			prevClip.addEventListener(MouseEvent.CLICK, doPrevious);
			pageBox.addChild(prevClip);
			//add the dividing line
			var headerLine:Sprite = new Sprite();
			headerLine.graphics.beginFill(0x999999, 1);
			headerLine.graphics.drawRect(0, 0, panelWidth - 5, 1);
			headerLine.graphics.endFill();
			headerLine.x = pageBox.x + 2.5;
			headerLine.y = pageBox.y + pageBox.height + 5;
			container.addChild(headerLine);
			/*var pageInput:TextField = new TextField();
			pageInput.type = TextFieldType.INPUT;
			pageInput.x = curWidth + spacer;
			pageInput.border = true;
			pageInput.width = 10;
			pageInput.height = 19;
			pageInput.name = 'pageInput';
			pageInput.addEventListener(KeyboardEvent.KEY_DOWN, checkKey);
			pageBox.addChild(pageInput);
			*/
			//curWidth += pageInput.width + spacer;
			//add the next button display
			var nextClip = new TextClip('next >');
			nextClip.x = curWidth + spacer;
			nextClip.name = 'nextClip';
			curWidth += nextClip.width + spacer;
			nextClip.addEventListener(MouseEvent.CLICK, doNext);
			pageBox.addChild(nextClip);
			searchLink = new SearchLink();
			//add the search link
			if(fromSearch == false){
				searchLink.x = panelWidth - searchLink.width;
				searchLink.setText('search');
				searchLink.buttonMode = true;
				searchLink.mouseChildren = false;
				searchLink.addEventListener(MouseEvent.CLICK, showSearch);
				searchLink.addEventListener(MouseEvent.ROLL_OVER, searchLink.onMouseOver);
				searchLink.addEventListener(MouseEvent.ROLL_OUT, searchLink.onMouseOut);
				pageBox.addChild(searchLink);
			}
			else{
				searchLink.x = panelWidth - searchLink.width;
				searchLink.setText('return to project');
				searchLink.buttonMode = true;
				searchLink.mouseChildren = false;
				searchLink.addEventListener(MouseEvent.CLICK, returnFromSearch);
				searchLink.addEventListener(MouseEvent.ROLL_OVER, searchLink.onMouseOver);
				searchLink.addEventListener(MouseEvent.ROLL_OUT, searchLink.onMouseOut);
				pageBox.addChild(searchLink);
			}
		}
		
		private function showSearch(e:MouseEvent){
			//draw the search box
			var searchHeight = 100;
			var searchBox:MovieClip = new MovieClip();
			searchBox.name = 'searchBox'
			searchBox.graphics.beginFill(0xFFFFFF, 1);
			searchBox.graphics.drawRect(0, 0, 300, searchHeight);
			searchBox.graphics.endFill();
			searchBox.y = 50;
			searchBox.x = 50;
			searchBox.graphics.beginFill(0x000000, 1);
			searchBox.graphics.drawRect(0, 0, 300, 1);
			searchBox.graphics.drawRect(0, searchHeight, 300, 1);
			searchBox.graphics.drawRect(0, 0, 1, searchHeight);
			searchBox.graphics.drawRect(300, 0, 1, searchHeight);
			searchBox.graphics.endFill();
			container.addChild(searchBox);
			var idText = new TextClip("Enter your search terms:");
			idText.x = 5;
			idText.y = 10;
			searchBox.addChild(idText);
			var searchText:TextField = new TextField();
			searchText.type = TextFieldType.INPUT;
			searchText.border = true;
			searchText.borderColor = 0x000000;
			searchText.width = 150;
			searchText.height = 20;
			searchText.x = 15;
			searchText.y = idText.y + idText.height + 10;
			searchText.name = 'searchText';
			searchText.addEventListener(KeyboardEvent.KEY_DOWN, checkKey);
			searchBox.addChild(searchText);
			var goText:TextClip = new TextClip('Go');
			goText.x = searchText.x + searchText.width + 10;
			goText.y = searchText.y;
			goText.addEventListener(MouseEvent.CLICK, startSearch);
			searchBox.addChild(goText);
			var closeText:TextClip = new TextClip('close');
			closeText.x = searchBox.width - closeText.width - 5;
			closeText.y = searchBox.height - closeText.height - 5;
			closeText.addEventListener(MouseEvent.CLICK, closeSearch);
			searchBox.addChild(closeText);
		}
		
		private function startSearch(e:MouseEvent){
			//do the search
			this.doSearch();
		}
		
		private function doSearch(){
			searchString = TextField(MovieClip(container.getChildByName('searchBox')).getChildByName('searchText')).text;
			if(searchString != ""){
				clearChildren('sequenceBox');
				clearChildren('pageBox');
				//add the loading movie, with the colors Duke, the colors!
				var ldMovie:LoadingMovie = new LoadingMovie(15, "Searching");
				ldMovie.x = (panelWidth - ldMovie.width)/2;
				ldMovie.y = 100;
				ldMovie.name = 'ldMovie';
				ldMovie.startAnimation();
				container.addChild(ldMovie);
				responder = new Responder(handleCount);
				connection = new NetConnection();
				connection.connect(MovieClip(root).gateway);
				if(curType == 'project'){
					query = "SELECT seq_id, project_id, description, type FROM " + curID + "_sequences WHERE description LIKE '%" + searchString + "%'";
				}
				else{
					query = "SELECT seq_id, project_id, description FROM " + curID + "_sequences WHERE description LIKE '%" + searchString + "%' and file_id='" + curFileID + "'";
				}
				connection.call("DBI.execSelect", responder, query);
				var searchBox:MovieClip = MovieClip(container.getChildByName('searchBox'));
				container.removeChild(searchBox);
				var tmp:TextField = TextField(this.getChildByName('description_text'));
				tmp.text = "Search Results For: " + searchString;;
				fromSearch = true;
				baseQuery = query;
			}
			else{
				this.closeSearch(new MouseEvent(MouseEvent.CLICK));
			}
		}
		
		private function returnFromSearch(e:MouseEvent):void{
			fromSearch = false;
			//load the project from where we left off
			if(curType == 'project'){
				this.loadSeqs(curID, curType, 1);
			}
			else{
				this.loadSeqs(curFileID, curType, 1);
			}
		}
		
		private function closeSearch(e:MouseEvent){
			var searchBox:MovieClip = MovieClip(container.getChildByName('searchBox'));
			container.removeChild(searchBox);
		}
		
		
		private function doPrevious(e:MouseEvent){
			this.loadPage(curID, curType, curPage - 1);
		}
		
		private function doNext(e:MouseEvent){
			this.loadPage(curID, curType, curPage + 1);
		}	
				
		private function checkKey(e:KeyboardEvent){
			if(e.charCode == 13){
				this.doSearch();
			}
		}
		
		private function clearChildren(clipName:String){
			lastY = 0;
			if(clipName == 'sequenceBox'){
				while(sequenceBox.numChildren > 0){
					/*var tmp = sequenceBox.getChildAt(0);
					if(flash.utils.getQualifiedClassName(tmp) == 'sequencepanel::Sequence'){
						container.removeEventListener(ExpandEvent.EXPAND_EVENT, tmp.onExpand);
					}*/
					sequenceBox.removeChildAt(0);
				}
			}
			else{
				while(MovieClip(container.getChildByName(clipName)).numChildren > 0){
					MovieClip(container.getChildByName(clipName)).removeChildAt(0);
				}
			}
		}
		
		private function addDescriptionText(str:String, type:String){
			//check if there already exists description text
			if(this.getChildByName('description_text')){
				this.removeChild(this.getChildByName('description_text'));
			}
			var descFormat:TextFormat = new TextFormat();
			var descFont:Font = new Font1();
			descFormat.font = descFont.fontName;
			descFormat.size = 14;
			descFormat.color = 0x025976;
			var descriptionText:TextField = new TextField();
			descriptionText.defaultTextFormat = descFormat;
			descriptionText.embedFonts = true;
			descriptionText.antiAliasType = AntiAliasType.ADVANCED;
			descriptionText.x = 0;
			descriptionText.text = str;
			descriptionText.name = 'description_text';
			descriptionText.autoSize = TextFieldAutoSize.LEFT;
			descriptionText.y = container.y - descriptionText.height - 5;
			this.addChild(descriptionText);				
		}
		
	}



}