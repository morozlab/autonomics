// JavaScript Document

function toggleDetail(target){
	var el = document.getElementById(target + "_display");
	var header = document.getElementById(target + "_header");
	if(el.style.display == 'none'){
		el.style.display = 'block';
	}
	else{
		el.style.display = 'none';
	}
}

function displayAnnotationAlignment(annotationID, sb, source, description){
	var request = createAjaxRequest();
	var alignmentContainer = document.getElementById(annotationID + "_display");
	var alignmentLink = document.getElementById(annotationID + "_link");
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			alignmentContainer.innerHTML = request.responseText;
			alignmentLink.innerHTML = "<a href='#' onclick=\"hideAnnotationAlignment(this.id, '" + sb + "', '" + source + "', '" + description + "');return(false);\" TITLE='View Annotation Alignment' id='" + annotationID + "'>" + description + "</a>";
			alignmentContainer.style.display = 'block';
		}
	}
	request.open("GET", "ajax/seq_detail/getAlignment.php?sb=" + sb + "&annotationID=" + annotationID + "&source=" + source, true);
	request.send(null);
}

function hideAnnotationAlignment(annotationID, sb, source, description){
	var alignmentLink = document.getElementById(annotationID + "_link");
	alignmentLink.innerHTML = "<a href='#' onclick=\"displayAnnotationAlignment(this.id, '" + sb + "', '" + source + "', '" + description + "');return(false);\" TITLE='View Annotation Alignment' id='" + annotationID + "'>" + description + "</a>";
	var alignmentContainer = document.getElementById(annotationID + "_display");
	alignmentContainer.style.display = 'none';
}
