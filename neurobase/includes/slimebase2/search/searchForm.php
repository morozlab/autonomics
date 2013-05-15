<form name='search' id='searchForm' method='get' action='viewSeqs.php'>
            	search <input type='hidden' name='ref' value='browse' /> 
                <select name='projectID' style='width: 150px;'>
                <option value='-1' selected="selected">Select Database</option>
				<?php
				if(!(isset($projectID))){
					$projectID = -1;
                }
				if($universalAccess == false){
					if($groupMembership != ""){
						$query = "SELECT projectID, project_name FROM project_directory, project_group_membership, user_group_membership WHERE user_id='" . $_SESSION['user_id'] . "' AND user_group_membership.group_id = project_group_membership.group_id AND project_directory.projectID = project_group_membership.project_id";
					}
					else{
						$query = "SELECT projectID, project_name, user_name FROM project_directory, permissions, users WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) AND project_directory.projectID = permissions.project_id AND permissions.user_id = '" . $_SESSION['user_id'] . "' AND users.user_id = permissions.user_id ORDER BY project_name";
					}
				}
				else{
					$query = "SELECT projectID, project_name FROM project_directory WHERE ((num_AA_seqs != 0) OR (num_NT_seqs != 0)) ORDER BY project_name";
				}
					$result = mysql_query($query);
					while($row = mysql_fetch_array($result)){
						$toEcho = "<option value='" . $row['projectID'] . "'";
						if($row['projectID'] == $projectID){
							$toEcho .= " selected='selected'";
						}
						$toEcho .= ">" . $row['project_name'] . "</option>";
						echo($toEcho);
					}
				?>
                </select>
                for <input type='text' size='20' name='query' /> <input type='submit' name='submit' value='Go' />
</form>