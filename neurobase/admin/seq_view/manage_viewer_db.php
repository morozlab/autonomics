<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head><?php
include("../../includes/checkSession.php");
include("../../includes/dbConnect.php");
include("src/recursiveDisplay.php");
mysql_select_db("moroz_lab");

//check to see if user submitted files
if(isset($_GET['mode']) && isset($_POST['uploadingFile'])){
	//if size == 0, no file was selected for upload
	if($_FILES['uploadedFile']['size'] == 0){
		echo($_FILES['uploadedFile']['size']);
		echo("<script language = 'javascript'>alert('Please select a file to upload!')</script>");
	}
	else{
		$projID = $_POST['projID'];
		$query = "SELECT path FROM project_directory WHERE projectID = '" . $projID . "'";
		$result = mysql_query($query);
		while($row = mysql_fetch_array($result)){
			$path = $row['path'];
		}
		$path = $path . "\\" . basename($_FILES['uploadedFile']['name']);
		if(move_uploaded_file($_FILES['uploadedFile']['tmp_name'], $path))
		{
			//insert entry into file list for that project
			$query = "INSERT INTO project_files SET type='" . $_POST['type'] . "_sequence', projectID ='" . $projID . "', file_size ='" .$_FILES['uploadedFile']['size'] . "', file_name = '" . basename($_FILES['uploadedFile']['name']) . "'";
			$result = mysql_query($query); 
			//get the fileID created for this file
			$fileID = mysql_insert_id();
			//fork perl script to insert all of the sequences into the database
			//$exe ='perl C:/wamp/www/admin/seq_view/src/scripts/insertSeqs.pl';
			//$args = '"' . $path . '" "' . $projID . '" "' . $fileID . '" "' . $_POST['type']; 
			//pclose(popen("start \"bla\" \"" . $exe . "\" " . escapeshellarg($args), "r"));
			shell_exec('perl C:/wamp/scripts/insertSeqs.pl "' . $path . '" "' . $projID . '" "' . $fileID . '" "' . $_POST['type'] . '" > dummy.txt');
			echo('perl C:/wamp/scripts/insertSeqs.pl "' . $path . '" ' . $projID . " " . $fileID . " " . $_POST['type'] . " > dummy.txt");		
		}
		else {
		echo "<script language='javascript'>alert('There was an error uploading your file. Please try again.');</script>";
		}
	}
}
?>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Manage SeqView DB</title>

<script type="text/javascript">
<!--
function MM_preloadImages() { //v3.0
  var d=document; if(d.images){ if(!d.MM_p) d.MM_p=new Array();
    var i,j=d.MM_p.length,a=MM_preloadImages.arguments; for(i=0; i<a.length; i++)
    if (a[i].indexOf("#")!=0){ d.MM_p[j]=new Image; d.MM_p[j++].src=a[i];}}
}
//-->
</script>
<link href="../../css/main.css" rel="stylesheet" type="text/css" />
</head>

<body class="mainPage" onload="MM_preloadImages('../images/layout/rollovers/home_over.jpg','../images/layout/rollovers/background_over.jpg','../images/layout/rollovers/pre_over.jpg','../images/layout/rollovers/post_over.jpg','../images/layout/rollovers/regulatory_over.jpg','../images/layout/rollovers/proposal_over.jpg','../images/layout/rollovers/moroz_lab_member_roll.jpg')">
<div id="container">
  <div id="header"></div>
   <?php include('../includes/admin_nav.php');?>
  <div id="content">
  <div id="leftContent">
        <p>
        <table width="390">
          <tr>
            <td width="148"><a href="manage_viewer_db.php?mode=add_remove">add/remove projects</a></td>
            <td width="113"><a href="manage_viewer_db.php?mode=link_files">add files to project</a></td>
            <td width="69">&nbsp;</td>
            <td width="40">&nbsp;</td>
          </tr>
        </table>
        </p>
        <h1>Current Project Structure</h1>
        <ul>
        	<li>Root</li>
		<?php recursiveDisplay(1); ?>
        </ul>
    </div>
<div id="rightContent">
   	    <?php 
		if(isset($_GET['mode'])){
			$mode = $_GET['mode'];
		}
		else{
			$mode = "add_remove";
		}
			if(($mode == 'add_remove') || (!(isset($mode)))){
		?>
        <h1>Add a New Project</h1>
   	    <form name="insert_project" method="POST" action="src/doInsertProject.php">
       	  <br />
       	  Project Name (unique)
		<input type="text" id="projName" name="projName" size="25" />
        	<br />
            <br />
        Insert as subdirectory of:<br />
        <select style="font-family: Verdana, Arial, Helvetica, sans-serif; font-size: 8pt;" name="projParent" id="projParent">
		<option value="1">root</option>
		<?php 
			$query = "SELECT * FROM project_directory WHERE project_name != 'root' ORDER BY project_name";
			$result = mysql_query($query);
			while($row = mysql_fetch_array($result)){
				$projectID = $row['projectID'];
				$projectName = $row['project_name'];
				$query2 = "SELECT project_name FROM project_directory WHERE projectID = '" . $row['child_of'] . "'";
				$result2 = mysql_query($query2);
				while($row2 = mysql_fetch_array($result2)){
					$toEcho = "<option value='" . $row['projectID'] . "'>" . $row['project_name'] .  " in " . $row2['project_name'] . "</option>";
					echo($toEcho);
				}
			}
		?>
        </select>
        <br />
        <br />
        <input type='submit' name='submit' id='submit' value='Go' />
      </form>
   	    <br />
   	    <br />
        <h1>Delete a Project</h1>
      <form name="delProj" id="delProj" method="POST" action="src/doDeleteProject.php">
        <br />
        Select a project to delete:<br />
        <select name="projID" id="projID">
        <?php 
			$query = "SELECT project_name, projectID FROM project_directory WHERE project_name != 'root' ORDER BY project_name";
			$result = mysql_query($query);
			while($row = mysql_fetch_array($result)){
				$projectName = $row['project_name'];
				$projectID = $row['projectID'];
				$toEcho = "<option value='" . $projectID . "'>" . $projectName . "</option>";
				echo($toEcho);
			}
		?>
        </select><br />
        <input type="submit" id="submit" name="submit" value="Go" />
      </form>
      <?php }
	  else{
	  ?>
      <h1>Add Sequences to a Project</h1>
      <br />
      <form method="post" enctype="multipart/form-data" action="manage_viewer_db.php?mode=link_files" name="linkFile">
      <p>Select a Project<br />
      <select name="projID" id="projID">
      <?php
	  	$query = "SELECT project_name, projectID, child_of FROM project_directory WHERE projectID != '1' ORDER BY project_name";
		$result = mysql_query($query);
		while($row = mysql_fetch_array($result)){
			$query2 = "SELECT project_name FROM project_directory WHERE projectID = '" . $row['child_of'] . "'";
			$result2 = mysql_query($query2);
			while($row2 = mysql_fetch_array($result2)){
				$toEcho = "<option value='" . $row['projectID'] . "'>" . $row['project_name'] .  " in " . $row2['project_name'] . "</option>";
				echo($toEcho);
			}
		}			
	  ?>
      </select><br /><br />
      Browse for Sequence File<br />
      <input type="file" name="uploadedFile" />
      <br /><br />
      Seuqneces in this file are: <select name='type' id='type'><option value='AA'>Protein</option><option value='NT'>Nucleotide</option></select>
      <input type="hidden" name="uploadingFile" value="yes" />
      <br />
      <input type="submit" value="Go" name="submit" />
      </p>
      </form>
      <?php } ?>   
   	</div>
  	</div>
	<div id="footer"></div>
</div>
</body>
</html>