<img style='float:left; margin-right: 20px;' src="images/slimeicon.jpg" alt="SlimeBase" />
        <div class='navElement'><a href='../index.php'>Home</a></div>
        <div class='navElement'>
        <?php
			if(isset($groupMembership) && ($groupMembership != "")){
				if($groupMembership == 2){
					?> <a href='nar.php'>Transcriptomes</a> <?php	
				}
				else if($groupMembership == 1){
					?> <a href='chpnas.php'>Transcriptomes</a><?php	
				}
			}
			else{ ?>
        		<a href='browse.php'>Transcriptomes</a>
        	<?php } ?>
        </div>
        <div class='navElement'><a href='blast.php?view=blast'>BLAST</a></div>

