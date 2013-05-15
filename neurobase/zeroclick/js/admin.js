var hidden = ['data-selection', 'input-data-select', 'project-tasks-header', 'project-tasks-checks', 'new-project-name', 'valid-project-name', 'configuration-submit'];
var valid_project_name = false;
var current_pn_mapping = "";

window.onload = function(){
	
	var sbmt = document.getElementById('configuration-submit');
	sbmt.onclick = submit_configuration;
	var data_select = document.getElementById('data-source-select');
	data_select.onchange = update_source_input;
	var proj_name = document.getElementById('input-project-name');
	proj_name.onblur = check_proj_name;
	document.getElementById('input-data-select').onchange = run_selection_changed;
	for(var i = 0; i < hidden.length; i++){
		set_class(hidden[i], "");
		append_class(hidden[i], "hidden");	
	}
}

function check_proj_name(){
	
	var proj_name = document.getElementById('input-project-name').value;
	var req = create_ajax_request();
	var source_id = document.getElementById('data-source-select').value;
	var valid = document.getElementById('valid-project-name');
	if(proj_name == ''){
		valid.innerHTML = "Project name cannot be empty.";
		valid_project_name = false;
		return;
	}
	req.onreadystatechange = function(){
		if(req.readyState == 4){
	//		alert(req.responseText)
			var ret_val = parseFloat(req.responseText);
			if(ret_val != 0 && proj_name != current_pn_mapping){
				valid_project_name = false;
				valid.innerHTML = "Project name taken.";
			}
			else{
				valid.innerHTML = "";
				valid_project_name = true;	
			}
		}
	}
	
	valid_project_name = false;
	req.open("GET", "/zeroclick/ajax/project.php?method=check_name&name=" + proj_name + "&source_id=" + source_id, true);
	req.send(null);
}

function init_proj_config(){
	
	var source_select = document.getElementById('data-source-select');
	var data_select = document.getElementById('input-data-select');
	var req = create_ajax_request();
	var run_name = data_select.options[data_select.selectedIndex].value;
	var source_id = source_select.options[source_select.selectedIndex].value;
	
	req.onreadystatechange = function(){
		
		if(req.readyState == 4){		
			var ret = JSON.parse(req.responseText);
			var checks = document.getElementById('project-tasks-checks');
			checks.innerHTML = ret['config'];
			var proj_id = document.getElementById('project-id');
			proj_id.value = ret['proj_id'];
			var proj_name = document.getElementById('input-project-name');
			proj_name.value = ret['proj_name'];
			current_pn_mapping = ret['proj_name'];
			valid_project_name = true;
		}
		
	}
	
	req.open("GET", "/zeroclick/ajax/project.php?method=init_config&run_name=" + run_name + "&source_id=" + source_id, true);
	req.send(null);
	
	
}

function run_selection_changed(){
	var statusBx = document.getElementById('proj-submit-progress').innerHTML = "";
	init_proj_config();
	
}

function submit_configuration(){
	
	//check if the user entered a project name
	var proj_el = document.getElementById('input-project-name');
	var proj_name = proj_el.value;
	if(proj_name == ''){
		proj_el.focus();
		var valid = document.getElementById('valid-project-name');
		valid.innerHTML = "Project name cannot be blank.";
		return;	
	}
	
	if(!valid_project_name){
		proj_el.focus()
		return;	
	}
	
	var req = create_ajax_request();
	var statusBx = document.getElementById('proj-submit-progress');
	
	req.onreadystatechange = function(){
		if(req.readyState == 4){
			var ret = parseFloat(req.responseText);
			if(ret == 0){
				statusBx.innerHTML = "Confguration successfully submitted!";
			}
			else{
				statusBx.innerHTML = "Error during configuration submission!";	
			}
		}
	}
	
	var form = document.getElementById('new-project-details');
	var getStr = "&" + form.elements[0].name + "=" + form.elements[0].value;
	
	for(i = 1; i < form.elements.length; i++){
		var el = form.elements[i];
		var val;
		if(el.type == 'select-one'){
			val = el.options[el.selectedIndex].value;
		}
		else if(el.type == 'checkbox'){
			val = el.checked;	
		}
		else{
			val = el.value;
		}	
		getStr += "&" + el.name + "=" + val;
			
	}
	req.open("GET", "/zeroclick/ajax/project.php?method=write_config" + getStr, true);
	req.send(null);
	
	statusBx.innerHTML = "Submitting configuration...";
	
}

function update_source_input(){
	
	var source_select = document.getElementById('data-source-select');
	var selected = source_select.options[source_select.selectedIndex].value;
	var data_div = document.getElementById('data-selection');
	var data_selection = document.getElementById('input-data-select');
	
	var req = create_ajax_request();
	
	req.onreadystatechange = function(){
		
		if(req.readyState == 4){
			var returnArray = JSON.parse(req.responseText);
			if(returnArray.type == 'ion_torrent' || returnArray.type == 'ftp' || returnArray.type == 'ion_proton' || returnArray.type == 'miseq'){
				data_selection.innerHTML = returnArray.option;
				if(returnArray.num_opts > 0){
					init_proj_config();
					for(var i = 0; i < hidden.length; i++){
						remove_class(hidden[i], "hidden");	
					}
				}
				else{
					remove_class('data-selection', 'hidden');
					remove_class('input-data-select', 'hidden');
				}
			}
		}
			
	}
	
	req.open("GET", "/zeroclick/ajax/data_sources.php?method=getInput&projectID=-5&src=" + selected, true);
	req.send(null);
		
}