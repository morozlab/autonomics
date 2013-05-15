<?php

include_once('include/std_include.php');
include_once($path . '/zeroclick/ajax/include/str.php');
include_once('../config/path.php');
include_once($path . '/zeroclick/ajax/include/zeroclick_database.php');

$method = $_GET['method'];

//$job_to_div = array('blast_swissprot' => 'blast-swissprot-task', 'blast_nr' => 'blast-nr-task', 'pfam' => 'pfam-annotation', 'go' => 'go-annotation', 'kegg' => 'kegg-annotation', 'upload' => 'upload-task');
$check_labels = array('assemble' => 'Assemble Transcripts', 'blast_swissprot' => 'SwissProt Annotation', 'blast_nr' => 'NR Annotation', 'pfam' => 'Pfam Domains', 'go' => 'Gene Ontology Terms', 'kegg' => 'KEGG Pathways', 'upload' => 'Upload to Neurobase', "quantification" => "Transcript Quantification", 'adapter_trim'=>'Adapter Trimming', 'quality_trim' => 'Quality Trimming', 'read_normalization' => 'Normalize Reads');
//array that maps input fields from the form to job_names in the zero-click system
//$job_types = array('assemble', 'blast_swissprot', 'blast_nr', 'pfam', 'go', 'kegg', 'upload');

if($method == 'write_config'){
	
	
	/*
	       <div id='new-project-name'>Name this project: <input name='input-project-name' id='input-project-name' type='text' /></div>
       <div id='project-tasks-header'>Select Tasks to Run:</div>
       <div id='project-tasks-checks'>
        	<div class='task-check'>BLAST against SwissProt <input name='blast-swissprot-task' type='checkbox' checked='checked' /><br /></div>
            <div class='task-check'>BLAST against NR <input name='blast-nr-task' type='checkbox' checked='checked' /><br /></div>
            <div class='task-check'>Pfam annotation <input name='pfam-annotation' type='checkbox' checked='checked' /><br /></div>
            <div class='task-check'>Gene Ontology Terms <input name='go-annotation' type='checkbox' checked='checked' /><br /></div>
            <div class='task-check'>KEGG <input name='kegg-annotation' type='checkbox' checked='checked' /><br /></div>
            <div class='task-check'>Upload to Neurobase <input name='upload-task' type='checkbox' checked='checked'  /><br /></div>
		*/
		
	//insert project name into pn_mapping and runname_to_projectname
	$project_name = preg_replace("[^A-Za-z0-9]", "_", $_GET['input-project-name']);
	$source_id = sanitize_number($_GET['data-source-select']);
	$pid = sanitize_number($_GET['project-id']);
	if($source_id == null){
		return;	
	}
	$run_name = $_GET['input-data-select'];
	
	$sth = $dbh->prepare("UPDATE pn_mapping SET project_name = ? WHERE project_id= ?");
	$sth->execute(array($project_name, $pid));		
	
	$sth = null;
	//gather the configuration bits to add to the db
	$select_config = $dbh->prepare("SELECT job_type FROM default_configuration");
	$select_config->execute(array());
	while($row = $select_config->fetch()){
		$job_type = $row['job_type'];
		if(isset($_GET[$job_type . '-checkbox'])){
			$checked = sanitize_checkbox($_GET[$job_type . '-checkbox']);
			$config = "-";
			if($checked){
				$config = "+";	
			}
			//check if this job has a configuration
			if(!job_exists($pid, $job_type)){//need to create a job id
				$sth = $dbh->prepare("INSERT INTO jn_mapping (project_id, job_type) VALUES (?, ?)");
				$sth->execute(array($pid, $job_type));
				$jid = $dbh->lastInsertId();
			}
			if(!job_configured($pid, $job_type)){
				$sth = $dbh->prepare("INSERT INTO configuration (project_id, job_type, code) VALUES (?, ?, ?)");
				$sth->execute(array($pid, $job_type, $config));
			}
			else{
				//update this job's configuration
				$sth = $dbh->prepare("UPDATE configuration SET code = ? WHERE project_id = ? AND job_type = ?");
				$sth->execute(array($config, $pid, $job_type));					
			}
			
		}
	}
	///update this project as configured
	$update = $dbh->prepare("UPDATE runname_to_pid SET configured='Y' WHERE project_id = ?");
	$update->execute(array($pid));
	echo(0);
}

else if($method == 'check_name'){
	
	$proj_name = trim($_GET['name']);
	//check list of mapped project names to see if there are any clashes
	$sth = $dbh->prepare("SELECT * FROM pn_mapping WHERE project_name = ?");
	$sth->execute(array($proj_name));
	if($sth->rowCount()){
		echo(-1);
		return;	
	}
	
	//if no clash, return status code 0
	echo(0);
		
}

else if($method == 'init_config'){
	
	$source_id = sanitize_number($_GET['source_id']) ;
	$run_name = $_GET['run_name'];
	
	//need to gather the information to display to the user and return it in an associative array
	//ret['proj_name'] = current project name, ret['config'] = configuration options
	$ret = array();
	
	//first, select the project id based on the run name and source id
	$sth = $dbh->prepare("SELECT pn_mapping.project_id, project_name FROM runname_to_pid JOIN pn_mapping ON runname_to_pid.project_id = pn_mapping.project_id WHERE run_name = ? and source_id = ?");
	$sth->execute(array($run_name, $source_id));
	$row = $sth->fetch();
	$ret['proj_name'] = $row['project_name'];
	$ret['proj_id'] = $row['project_id'];
	$ret['config'] = "";
	//select the configuration for this project
	$sth = $dbh->prepare("SELECT * FROM configuration WHERE project_id=?");
	$sth->execute(array($ret['proj_id']));
	$configs = array();
	while($row = $sth->fetch()){
		$configs[$row['job_type']] = $row['code'];
	}
	$sth = $dbh->prepare("SELECT * FROM default_configuration");
	$sth->execute();
	while($row = $sth->fetch()){
		$job_type = $row['job_type'];
		$checkbx = "<div class='task-check'>" . $check_labels[$job_type] . "<input type='checkbox' name='" . $job_type . "-checkbox' id='"  . $job_type . "-checkbox' ";
		$checkbx_end = " /></div>";
		$checked = '';
		if(isset($configs[$row['job_type']])){
			if($configs[$row['job_type']] != "-"){
				$checked = "checked='CHECKED'";	
			}
		}
		$ret['config'] .= $checkbx . $checked . $checkbx_end;
	}
	echo(json_encode($ret));
	
}

?>