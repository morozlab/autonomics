<?php 
session_start();
include('../includes/dbConnect.php');
include('../includes/restrictAccess.php');
include('../includes/slimebase2/projects/checkProjectResources.php');
mysql_select_db('moroz_lab');
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<link rel="icon" type="image/png" href="../favicon.png" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>NeuroBase: A Comparative Neurogenomics Database</title>
<link href="../css/neurobase_basic.css" rel="stylesheet" type="text/css" />
<link href="../css/browse.css" rel="stylesheet" type="text/css" />
<link href='../css/sequence_tab.css' rel='stylesheet' type='text/css'/>
<link href='../css/pfam_tab.css' rel='stylesheet' type='text/css' />

<script language='javascript' type='text/javascript' src="javascript/shadedborder/shadedborder.js"></script>
<script type='text/javascript' src='javascript/pfam/pfamUtils.js'></script>
<script language='javascript' type='text/javascript' src='javascript/mootools/mootools.js'></script>
<script language='javascript' type='text/javascript' src='javascript/projects/projectUtils.js'></script>
<script language='javascript' type='text/javascript' src="javascript/seqUtils.js"></script>

 <!--[if IE]>
      <script type="text/javascript" src="http://pfam.sanger.ac.uk/static/javascripts/excanvas.js"></script>
      <script type="text/javascript" src="http://getfirebug.com/releases/lite/1.2/firebug-lite-compressed.js"></script>
 <![endif]-->
 

<script language='javascript' type='text/javascript'>
	var detailExpanded = false;
	var currentSelected = "";
	sfHover = function() {
		var sfEls = document.getElementById("projectOptions").getElementsByTagName("LI");
		for (var i=0; i<sfEls.length; i++) {
			sfEls[i].onmouseover=function() {
				this.className+=" sfhover";
			}
			sfEls[i].onmouseout=function() {
				this.className=this.className.replace(new RegExp(" sfhover\\b"), "");
			}
		}
	}
	if (window.attachEvent) window.attachEvent("onload", sfHover);
</script>
</head>

<body>
<div id='neurobase-logo'><img src="images/neuro_base_logo.png" width="624" height="69" /></div>
        <div id='neurobase-nav-left'><img src="images/neuro_base_nav_left.png" /></div>
        <div id='neurobase-nav'><div class='first-nav'><a href='../index.php'>home</a></div><div class='nav-element'><?php
			if(isset($groupMembership) && ($groupMembership != "")){
				if($groupMembership == 2){
					?> <a href='nar.php'>transcriptomes</a> <?php	
				}
				else if($groupMembership == 1){
					?> <a href='chpnas.php'>transcriptomes</a><?php	
				}
			}
			else{ ?>
        		<a href='browse.php'>transcriptomes</a>
        	<?php } ?></div><div class='nav-element'><a href='blast.php?view=blast'>BLAST</a></div><div class='nav-element'><a href='viewSeqs.php?searchProj=-1'>search</a></div></div>
<div id='container'>
	<div id="content">
    <?php include('../includes/slimebase2/displayLogin.php'); ?>
    <div id='mainContent'>
    	<div id='selectContainer'>
            <div id='selectBoxMiddle'>
            	<div id='selectBox'>
                	<form id='projectSelect' name="projectSelect">
                    	<select id='projectSelectElement' name='projectSelectElement'>
                        	<option value='-1'>Select species & tissue</option>
                            <?php
	  $order_by = " ORDER BY display_group, display_order ";
								if($universalAccess == false){
			$query = "SELECT projectID, project_name FROM project_directory, permissions WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND project_directory.projectID = permissions.project_id AND permissions.user_id ='" . $_SESSION['user_id'] . "'" . $order_by;
				}
				else{
					$query = "SELECT projectID, project_name FROM project_directory WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0))" . $order_by;
				}
								$result = mysql_query($query);
								while($row = mysql_fetch_array($result)){
									print("<option value='" . $row['projectID'] . "'>" . $row['project_name'] . "</option>");
								} 
							?>
                        </select>
                    </form> 
                </div>
                <p style='margin-left:40px; color:#fff;'>This site is currently under construction.</p>
            </div>
        </div>
        <div id='projDetailContainer'>
        	<div id='projDetailContent'>
            	<div id='organism-description'></div>
            	<div id='projectImage'></div>
                <div id='projectStats'>
                    <div id='pSContent'><h1>Assembly Statistics</h1>
                    	<p id='pSContentParagraph'>
                        sequenced using 454<br />
                        last updated on 6/20/2010<br />
                        7658 contigs<br />
                        765 singlets<br />
                        303bp average contig length<br />
                        6000 annotated sequences<br />
                        </p>
                    </div>
                </div>
            </div>
        </div>
   </div>
   </div> 
</div>
<div id='main-information-content'>
<h1>Transcriptomes hosted on Neurobase were designed and constructed at the University of Florida</h1>
<div class='main-people-details'>Professor Leonid L. Moroz<br/>Princple Investigator<br/>9505 Oceanshore Blvd<br/>St Augustine, FL 32080<br/>904-461-4020<br/><a href='mailto:moroz@whitney.ufl.edu'>moroz@whitney.ufl.edu</a></div><div class='main-people-details'>Mathew Citarella<br/>Bioinformatics/Design<br/>9505 Oceanshore Blvd<br/>St Augustine, FL 32080<br/>904-461-4007<br/><a href='mailto:mathew.citarella@gmail.com'>mathew.citarella@gmail.com</a></div><div class='main-people-details'>Dr. Andrea Kohn<br/>Genomics/Molecular Biology<br/>9505 Oceanshore Blvd<br/>St Augustine, FL 32080<br/>904-461-4007<br/><a href='mailto:abkohn@msn.com'>abkohn@msn.com</a></div><div class='contact-message'><p>Please direct comments or questions about the database to <a href='mailto:morozgenomics@whitney.ufl.edu'>morozgenomics@whitney.ufl.edu</a>.</p></div>
</div>
<div id='extraContent'>
	<div id='extraContentNav'>
    </div>
    <div id='extraContentMain'>
    	
    </div>
</div>
<?php
if(isset($_SESSION['project_loaded'])){	
?>
<script type='text/javascript' language="javascript">
	var projSelect = document.getElementById('projectSelectElement');
	projSelect.selectedIndex = 0;
</script>
<?php
}
?>
</body>
</html>
