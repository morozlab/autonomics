// JavaScript Document

function addToCustomFasta(link_id){
	var idParts = link_id.split("_");
	var sb = idParts[0];
	var projectID = idParts[1];
	var request = createAjaxRequest();
	request.onreadystatechange = function(){
		if(request.readyState == 4){
			var response = request.responseText;
			//update the custom fasta counter with new number
			var el = document.getElementById("fastaCounter");
			el.innerHTML = "FASTA (" + response + ")";
			var container = document.getElementById(sb + "_customFastaLink");
			container.innerHTML = "<img border='0' id='" + sb + "_" + projectID + "_img' src='images/view_seqs/added_to_custom_fasta.jpg' title='Already Added to Custom Fasta' />";
		}
	}
	request.open("GET", "ajax/project_view/utils.php?method=updateCart&args=" + sb + "_" + projectID, true);
	request.send(null);
}