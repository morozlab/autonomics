<?php 
session_start();
include('../includes/dbConnect.php');
include('../amfphp/services/Util.php');
include('../includes/restrictAccess.php');
mysql_select_db('moroz_lab');
if(isset($_GET['sbid'])){
	$sbid = $_GET['sbid'];
}
else{
	$sbid = -1;
}

if(isset($_GET['projectID'])){
	$projectID = $_GET['projectID'];
}
else{
	$projectID = -1;
}
include('../includes/slimebase2/security/checkPermissions.php');
//get the project name
$query = "SELECT project_name FROM project_directory WHERE projectID ='" . $projectID . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$projectName = $row['project_name'];
//look up which project the sbid is FRO
$query = "SELECT NT_sequence, AA_sequence, sb_id, description, length FROM " . $projectID . "_sequences WHERE sb_id ='" . $sbid . "'";
$result = mysql_query($query);
$row = mysql_fetch_array($result);
$description = $row['description'];
$sb = $row['sb_id'];
$nt = $row['NT_sequence'];
$length = $row['length'];
if(is_null($row['AA_sequence'])){
	$aa = "translate";
}
else{
	$aa = $row['AA_sequence'];
}
if($row['AA_sequence'] == ""){
	$aa = "translate";
}

//get the top annotation
$query = "SELECT text FROM sorted_homology, annotation_db WHERE sorted_homology.sb_id = '" . $sb . "' AND sorted_homology.annotation_id = annotation_db.annotation_id AND sorted_homology.sort_id='1' AND sorted_homology.project_id='" . $projectID . "' LIMIT 0,1";
$result = mysql_query($query);
while($row = mysql_fetch_array($result)){
	$description = $row['text'];
} 


?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title><?php echo("Sequence Detail - sb|" . $sbid . "|"); ?></title>

<link rel="stylesheet" href="../css/pfam.css" type="text/css" />
<link rel="stylesheet" href="../css/prototip.css" type="text/css" />
<link href='../css/slimebase_basic.css' rel='stylesheet' type='text/css' />
<link href="../css/slimebaseSeqDetail.css" rel="stylesheet" type="text/css" />
</head>

<body>

<!-- load neurobase utility scripts -->
<script type='text/javascript' src="javascript/seqUtils.js"></script>
<script type='text/javascript' src='javascript/displayUtils.js'></script>
<!-- load the prototype library from google --> 
<script type="text/javascript" src="http://www.google.com/jsapi"></script> 
<script type="text/javascript">google.load("prototype", "1.7");</script> 
<!-- load pfam domain libraries -->
<script type='text/javascript' src='javascript/pfam/excanvas.js'></script> 
<script type='text/javascript' src='javascript/pfam/canvas.text.js?reimplement=true&amp;dontUseMoz=true'></script> 
<script type='text/javascript' src='http://pfam.sanger.ac.uk/static/javascripts/faces/optimer-bold-normal.js'></script> 
<script type='text/javascript' src='javascript/pfam/domain_graphics.js'></script> 
<script type='text/javascript' src='javascript/pfam/prettify.js'></script> 
<script type='text/javascript' src="javascript/pfam/prototip.js"></script>
<script type='text/javascript' src="javascript/pfam/styles.js"></script>
<script type='text/ecmascript' src='javascript/pfam/pfamRendering.js' />

 <!--[if IE]>
      <script type="text/javascript" src="http://pfam.sanger.ac.uk/static/javascripts/excanvas.js"></script>
      <script type="text/javascript" src="http://getfirebug.com/releases/lite/1.2/firebug-lite-compressed.js"></script>
 <![endif]--> 


<!-- set the projectID variable for use elsewhere -->
<script language='javascript' type='text/javascript'>
	var projectID = <?php echo($projectID); ?>;
</script>
<div id='container'>
	<div id="content">
	<div id='neuroBaseHeader'>
    	<img style='float:left; margin-right: 20px;' src="images/slimeicon.jpg" alt="SlimeBase" />
        <div class='navElement'><a href='../index.php'>Home</a></div>
        <div class='navElement'>
        <?php
			if(isset($groupMembership) && ($groupMembership != "")){
				if($groupMembership == 2){
					?> <a href='nar.php'>Projects</a> <?php	
				}
				else if($groupMembership == 1){
					?> <a href='chpnas.php'>Projects</a><?php	
				}
			}
			else{ ?>
        		<a href='browse.php'>Projects</a>
        	<?php } ?>
        </div>
        <div class='navElement'><a href='blast.php?view=blast'>BLAST</a></div>
    	<div id="search">
        	<?php include ('../includes/slimebase2/search/searchForm.php'); ?>
        </div>
    </div>
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='contentHeader'></div>
        <div id='meat'>
        	<div id='contentTop'>
            	<div id='topLeftHeader'></div>
                <div id='middleHeader'></div>
                <div id='topRightHeader'></div>
            </div>
            <div id='contentBody'>
           		<div class='detail'>
                <h1 id='nucleotideHeader'><?php 
					echo("Sequence Details" . " - sb|" . $sb . "|");				
				?></h1>
                <div class='detailElement' id='nulceotide_display'>
                <p><?php
				if($aa == "translate"){ 
					echo(">sb|" . $sb . "|<br/>");
					echo($nt);
				}
				else{
					echo(">" . $description . "<br/>");	
					echo($aa);
				}
			    ?></p></div>
                <div  class='detailElement' id='annotation_display'>
                    	<?php if(isset($description)){
							?>
                            <div class='annotation_element'>
                        	<h2>Public Database Homology</h2>
                            <table width=100% id='annotationTable'><tr class='annotationHeaderRow'><td class='sourceHeader'>source</td><td class='descriptionHeader'>description</td><td class='smallNumberHeader'>e-value</td></tr>
                    		<?php
                            	$query = "SELECT sorted_homology.annotation_id, evalue, name, text, sorted_homology.source FROM sorted_homology, annotation_db, annotation_sources WHERE annotation_db.annotation_id = sorted_homology.annotation_id AND annotation_db.source = annotation_sources.id AND sorted_homology.sb_id ='" . $sb . "' AND sorted_homology.sort_id='1' AND sorted_homology.project_id='" . $projectID . "'";
								$result = mysql_query($query);
								while($row = mysql_fetch_array($result)){
									$name = $row['name'];
									$evalue = $row['evalue'];
									$description = $row['text'];
									$id = $row['annotation_id'];
									$source = $row['source'];
									if($name == "sp"){
										$dbName = "SwissProt";
										$dbLink = "http://www.uniprot.org/uniprot/" . $id;
										$sourceLink = "<a href='http://www.uniprot.org/uniprot/" . $id . "' TITLE='View Uniprot Entry'><img border='0' src='images/view_seqs/Swiss-Prot-logo.png' /></a>";
									}
									else{
										$dbName = "NR";
										$dbLink = "http://www.ncbi.nlm.nih.gov/protein/" . $id;
										$sourceLink = "<a href='http://www.ncbi.nlm.nih.gov/protein/" . $id . "' TITLE='View NR Entry'><img border='0' src='images/view_seqs/white_ncbi.png' /></a>";
									}
									echo("<tr class='annotationDataRow'><td class='externalLinkCell'>");
									echo($sourceLink);
									echo("</td>");
									echo("<td class='annotationDescriptionCell'>");
									echo("<div class='annotationAlignmentLink' id='" . $id . "_link'><a href='#' onclick=\"displayAnnotationAlignment(this.id, '" . $sb . "', '" . $source . "', '" . $description . "');return(false);\" TITLE='View Annotation Alignment' id='" . $id . "'>" . $description . "</a></div><div class='annotationAlignmentDisplay' id='" . $id . "_display'></div>");
									echo("</td>");
									echo("<td class='evalCell'>");
									echo($evalue);
									echo("</td></tr>");
									
								}
							?>
                   			</table>
                        </div>
                        <?php }
						?>
                        <?php 
						//cehck for pfam annotations
						$query = "SELECT pfama.pfamA_acc, pfama.pfamA_id, description, evalue, start, stop, type FROM pfam_annotations, pfama WHERE project_id='" . $projectID . "' AND sb_id='" . $sbid . "' AND pfam_annotations.pfamA_acc = pfama.pfamA_acc";
						$result = mysql_query($query);
						if(mysql_num_rows($result)){
							?>
                            	<div class='annotation_element'>
                                	<h2>Pfam Domains</h2>
                                    <div id="pfam-detail-<?php echo($sbid); ?>-rendering"></div>
                                    <table width=100% id='pfamDomainTable' class='seqDetailTable'><tr class='annotationHeaderRow'><td class='sourceHeader'>domain id</td><td class='descriptionHeader'>description</td><td class='smallNumberHeader'>e-value</td><td class='smallNumberHeader'>start</td><td class='smallNumberHeader'>stop</td></tr>
                                    <script type='text/javascript' language="javascript">
                                    featureArray = new Array();
                                    </script>
                                    <?php
									$jsToExec = "";
									while($row = mysql_fetch_array($result)){
										$link = "<a href='http://pfam.sanger.ac.uk/family?acc=" . $row['pfamA_acc'] . "' target='_blank' TITLE='View Details For: " . $row['description'] . "'>" . $row['pfamA_acc'] . "</a>";
										echo("<tr class='annotationDataRow'><td class='externalLinkCell'>");
										echo($link);
										echo("</td>");
										echo("<td class='annotationDescriptionCell'>");
										echo($row['description']);
										echo("</td>");
										echo("<td class='evalCell'>");
										echo($evalue);
										echo("<td class='evalCell'>");
										echo($row['start'] . "</td>");
										echo("<td class='evalCell'>");
										echo($row['stop'] . "</td>");
										echo("</tr>");
										$jsToExec .= "featureArray.push(new pfamAnnotation('" . $row['pfamA_acc'] . "', '" . $row['start'] . "', '" . $row['stop'] . "', '" . $row['description'] . "', '" . $row['type'] . "', '" . $row['pfamA_id'] . "'));";
									}
									?>
                                    </table>
                                </div>
                            <?php	
							echo("<script type='text/javascript'>");
							echo($jsToExec);
							echo("displayPfamGraphic(featureArray, '" . $length . "', '" . $sbid . "');");
							echo("</script>");
						}
						
						//check for KEGG Pathways for this sequence
						$query = "SELECT * FROM kegg_annotations WHERE project_id='" . $projectID . "' AND sb_id='" . $sbid . "'";
						$result = mysql_query($query);
						if(mysql_num_rows($result) != 0){
						?>
                        
                        <div class='annotation_element'>
                        	<h2>KEGG Pathways</h2><table width=100% id='KeggPathwayTable' class='seqDetailTable'><tr class='annotationHeaderRow'><td class='sourceHeader'>pathway id</td><td class='descriptionHeader'>description</td><td></td><td class='smallNumberHeader'>e-value</td></tr>
                            <?php
							while($row = mysql_fetch_array($result)){
								$link = "<a href='http://www.kegg.jp/kegg-bin/search_pathway_text?map=map&keyword=" . trim($row['kegg_pathway_description']) . "&mode=1&viewImage=true' target='_blank' TITLE='View Details For: " . $row['kegg_pathway_description'] . "'>" . $row['kegg_id'] . "</a>";
								echo("<tr class='annotationDataRow'><td class='externalLinkCell'>");
								echo($link);
								echo("</td>");
								echo("<td class='annotationDescriptionCell'>");
								echo($row['kegg_pathway_description']);
								echo("</td>");
								echo("<td><a href='kegg.php?projectID=" . $projectID . "&pathwayID=" . $row['path_id'] . "'>view pathway diagram</a></td>");
								echo("<td class='evalCell'>");
								echo($row['evalue']);
								echo("</tr>");	
							}
						
							?>
                            </table>                            
                        </div>
                        <?php	
						}
						
						//see if there is a GO annotation for this sequence
							$query = "SELECT go_catalog.description, go_catalog.go_id, go_annotation_new.annotation_ev FROM go_annotation_new, go_catalog WHERE sb_id ='" . $sbid . "' AND go_catalog.go_id = go_annotation_new.go_id";
							$result = mysql_query($query);
							if(mysql_num_rows($result) != 0){
						?>                        
                        <div class='annotation_element'>
      						<h2>GO term(s)</h2>
                            
							<?php 
							while($row = mysql_fetch_array($result)){
							echo("<p><a href='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=" . $row['go_id'] . "' target='_blank' title='View This Term @ geneontology.org'>GO:" . $row['go_id'] . "</a> " . $row['description'] . "</p>"); 
							}
							?>
                        </div>
                        <?php } 
						?>
                 </div>
           	     </div>
                 <div class='detail'>
				 <?php if($aa == "translate"){
						echo("<div class='seqDetailHeaderBox'><h1 id='amino_acid_header'><a href='#' name='amino_acid' onclick=\"toggleDetail(this.name);return(false);\" target='_self' title='expand/collapse amino acid'>translate to amino acid</a></h1></div>");
						echo("<div class='detailElement' id='amino_acid_display' style='display:none;'></div><script language='javascript'>translateSequence('" . $nt . "');</script>");
					}
					?>
                  
                 </div>
                 	<?php //see if there is any expression data for this sequence
						$baseQuery = "SELECT type, description, filename FROM expression WHERE sb_id = '" . $sbid . "'";
						$result = mysql_query($baseQuery);
						if(mysql_num_rows($result)){
					?>
                 <div class='detail'>
                 	<div class='seqDetailHeaderBox'>
                    <h1 id='expression_header'><a href='#' name='expression' onclick="toggleDetail(this.name);return(false);" target='_self' title='expand/collapse expression data'>expression data</a></h1>
                    </div>
                    	<div class='detailElement' style='border:none;' id='expression_display'>
                        	<?php //check if we have in-situ expression data
								$query = $baseQuery . " AND type='IS'";
								$result = mysql_query($query);
								if(mysql_num_rows($result)){
							?>
                            <div class='annotation_element'>
                            	<h2>in-situ</h2>
                                <?php
								//iterate over the electrophys results and display each
								while($row = mysql_fetch_array($result)){
									$toEcho = "<div class='expressionImgDisplay'><img src='resources/expression/in-situ/" . $row['filename'] . "' alt='" . $row['description'] . "' /></div>";
									echo($toEcho);
								}
								?>
                            </div>
                            <?php
								}
								//check if we have electrophys expression data
								$query = $baseQuery . " AND type='EP'";
								$result = mysql_query($query);
								if(mysql_num_rows($result)){
							?>
                            <div class='annotation_element'>
                            	<h2>electrophysiology</h2>
                                <?php
								//iterate over the electrophys results and display each
								while($row = mysql_fetch_array($result)){
									$toEcho = "<div class='expressionImgDisplay'><img src='resources/expression/electrophysiology/" . $row['filename'] . "' alt='" . $row['description'] . "' /><div class='expressionCaption'><p>" . $row['description'] . "</p></div></div>";
									echo($toEcho);
								}
								?>
                            </div>
                            <?php }
							?>
                        </div>
                 </div>    	
                    <?php } ?>
            </div>
        	<div id='contentBottom'>
            	<div id='leftFooter'><img src="images/bottom_left.png" /></div>
                <div id='middleFooter'></div>
                <div id='rightFooter'><img src="images/bottom_right.png" /></div>
            </div>   
    	</div>
        <div id='sideNav'></div>
    </div>   
    </div>
</div>
</body>
</html>

<?php 
