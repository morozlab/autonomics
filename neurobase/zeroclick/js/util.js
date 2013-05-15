// JavaScript Document

function append_class(el_id, class_name){
	document.getElementById(el_id).className += " " + class_name;
}

function build_select_options(option_dict){
	ret = "";
	for(key in option_dict){
		ret += "<option value='" + key + "'>" + option_dict[key] + "</option>";
	}
	return ret;
}


function create_ajax_request(){
	var request;
	try{
		// Opera 8.0+, Firefox, Safari
		request = new XMLHttpRequest();
	} catch (e){
		// Internet Explorer Browsers
		try{
			request = new ActiveXObject("Msxml2.XMLHTTP");
		} catch (e) {
			try{
			request = new ActiveXObject("Microsoft.XMLHTTP");
			} catch (e){
				// Something went wrong
				alert("Your browser is out of date and does not support Ajax. Please update your browser to view GO Annotations");
				return false;
			}
		}
	}
	return request;
}

function remove_class(el_id, class_name){
	var re = new RegExp(" " + class_name + "\s*", "g");
	document.getElementById(el_id).className = document.getElementById(el_id).className.replace(re , '' );		
}

function set_class(el_id, class_name){
	document.getElementById(el_id).className = class_name;	
}
