<img style='float:left; margin-right: 20px;' src="images/slimeicon.jpg" alt="SlimeBase" />
        <div class='navElement'><a href='../index.php'>Home</a></div>
        <div class='navElement'><a href='browse.php'>Projects</a></div>
        <div class='navElement'><a href='blast.php?view=blast'>BLAST</a></div>
        <div class='navElement'>Quantify</div>
        <div class='navElement'><a href='fasta.php' id='fastaCounter'>FASTA
        <?php 
			$query = "SELECT COUNT(*) as C FROM sequence_cart WHERE session_id='" . session_id() . "'";
			$result = mysql_query($query);
			$num = 0;
			while($row = mysql_fetch_array($result)){
				$num = $row['C'];
			}
			echo("(" . $num . ")");
		?>
        </a>
        </div>
