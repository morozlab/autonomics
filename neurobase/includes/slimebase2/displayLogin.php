<div id='loginContainer'>
        	<div id='loginText'><?php
				$query = "SELECT name FROM users WHERE user_id ='" . $_SESSION['user_id'] . "'";
				$result = mysql_query($query);
				$row = mysql_fetch_array($result);
				$userName = $row['name'];
				if($_SESSION['user_id'] == 2){
					echo("viewing as <b>public</b> | <a href='login.php?ref=browse' title='log In'>log in</a>");
				}
				else{
					echo("logged in as <b>" . $userName . "</b> | <a href='logout.php' title='logout'>logout</a>");
				}
			?></div>
</div>
