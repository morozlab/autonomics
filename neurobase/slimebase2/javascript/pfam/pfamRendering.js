function displayPfamGraphic(featureArray, length, sbid){
	var intersections = new Array();
	var protein = '{"length":"' + length + '",';
	protein = decorateProtein(protein, featureArray);
	protein += "}";
	var sequence = eval('sequence = ' + protein);
	var pg = new PfamGraphic();
	pg.setParent("pfam-detail-" + sbid + "-rendering");
	try {
		var xscaling = 1500/length;
		var yscaling = xscaling;
		pg.setImageParams( {
          xscale: xscaling,
          yscale: 1.0
        } );
          pg.setSequence( sequence );
          pg.render();
    } catch ( e ) {
		alert(e);
          $("errors").update( e ).show();
          return;
    }
}

function decorateProtein(protein, features){
	protein += "regions:[";
	for(i = 0; i < features.length; i++){
		//generate the hex color for this feature
		var name = features[i].name;
		var rgbColors = colorHash(features[i]);
		var color = rgbToHex(rgbColors[0], rgbColors[1], rgbColors[2]);
		//add this feature to the protein
		
		//determine the edges for this domain
		var thisStart = parseInt(features[i].start);
		var thisEnd = parseInt(features[i].end);
		if(i != 0){
			var prevStart = parseInt(features[i - 1].start);
			var prevEnd = parseInt(features[i -1].end);
		}
		if(features[i + 1] != undefined){
			var nextStart = parseInt(features[i + 1].start);
			var nextEnd = parseInt(features[i + 1].end);	
		}
		if(features[i].type == 'Domain' || features[i].type == 'Family'){
			protein += '{"startStyle":';
			if(i != 0){
				if(prevEnd > thisStart){
					protein += '"jagged",';	
				}
				else{
					protein += '"curved",';	
				}
			}
			else{
				protein += '"curved",';	
			}
			if(features[i + 1] != undefined){
				if(nextStart < thisEnd){
					protein += '"endStyle":"jagged",';	
				}
				else{
					protein += '"endStyle":"curved",';	
				}	
			}
			else{
				protein += '"endStyle":"curved",';	
			}			
		}
		if(features[i].type == 'Motif' || features[i].type == 'Repeat'){
			protein += '{"startStyle":"straight", "endStyle":"straight",';	
		}
		protein += '"display":true,';
		//add the start and the end
		protein += '"start":"' + features[i].start + '", "end":"' + features[i].end + '",';
		protein += '"aliStart":"' + (parseInt(features[i].start) + 10) + '", "aliEnd":"' + (parseInt(features[i].end) - 10) + '",';
		//add the color
		protein += '"colour":"#' + color + '",';
		//add some metadata stuffs
		protein += '"metadata":{ "scoreName" : "e-value", "score" : "' + features[i].evalue + '", "description" : "' + features[i].description + '", "accession" : "'+ features[i].acc + '", "end" : "' + features[i].end + '", "database" : "pfam", "identifier" : "' + features[i].id + '", "type" : "' + features[i].type + '", "start" : "' + features[i].start + '"},';
		//add the text
		protein += '"text":"' + features[i].id.substr(0, 8) + '"';
		protein += '},';
	}
	protein += "]";
	return protein;
}


function pfamAnnotation(name, start, end, description, type, id, evalue, acc){
	this.name = name;
	this.start = start;
	this.end = end;
	this.description = description;
	this.type = type;
	this.id = id;
	this.evalue = evalue;
	this.acc = acc;
}

function rgbToHex(R,G,B) {return toHex(R)+toHex(G)+toHex(B)}

function toHex(n) {
 n = parseInt(n,10);
 if (isNaN(n)) return "00";
 n = Math.max(0,Math.min(n,255));
 return "0123456789ABCDEF".charAt((n-n%16)/16)
      + "0123456789ABCDEF".charAt(n%16);
}

function colorHash(feature){
	var name = feature.name;
	var value = 0;
	var colors = new Array();
	for (j=0;j < name.length; j++) {
    	value += parseInt(name.charCodeAt(j));
	}
	colors[0] = value % 256;
	value = 0; 
	for(j = 0; j < feature.description.length; j++){
		value += parseInt(feature.description.charCodeAt(j));
	}
	colors[1] = value % 256;
	colors[2] = (value * 3) % 256;
	return colors;
}
