<?php
session_start();
if(isset($_GET['msg'])){
	$msg = $_GET['msg'];
}
else{
	$msg = "none";
}
if(isset($_SERVER['HTTP_REFERER'])){
	$ref = $_SERVER['HTTP_REFERER'];
}
else{
	$ref = "http://74.252.103.104:8888/slimebase2/browse.php";
}	
if(isset($_GET['origref'])){
	$ref = $_GET['origref'];
}
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Log In to SlimeBase</title>
<link href="../css/slimebase_white.css" rel="stylesheet" type="text/css" />
<link href="../css/login.css" rel="stylesheet" type="text/css" />
</head>

<body>
<div id='container'>
	<div id="content">
	<div id='header'>
    	<?php include('../includes/slimebaseNav.php'); 
		
		?>
        
    <div id="search">
        </div>
    </div>
    <div id='mainContent' style='border:none; background-color:transparent;'>
    	<div id='errorMessage'>
        	<?php
				if($msg == 'permfail'){
					echo("We're sorry, but you don't have permission to view the selected project. Your login may have expired, or you may not have been granted access to the project. Please log in below.");
				}
				else if($msg == 'credentials'){
					echo("<script language='javascript'>alert('Incorrect login information! Please try again.');</script>");
				}
			?>
        </div>
    	<form class="membersArea" name="authForm" method="post" action="../lab_members/doLogin.php">
          <fieldset id="maForm">
   		  <legend>Log in</legend>
       	  <p><label for="userName">Username</label> <input type="text" name="userName" size="40"/></p>
			<p><label for="passWord">Password</label> <input type="password" name="passWord" size="40"/></p>
			<p class="submitButton"><input id="smButton" name="submit" type="submit" value="Log in" />
            <span style="margin-left: 1.5em; font-size: 8pt;"><a href="forgotten_password.php">forgot password</a> | <a href="forgotten_username.php">forgot username</a></span>			</p>
            </fieldset>
            <input type='hidden' name='ref' value='<?php echo($ref); ?>' />
            <?php if(isset($_GET['projectID'])){
				?>
				<input type='hidden' name='projectID' value='<?php echo($_GET['projectID']); ?>' />
            <?php 
			}
			?>
        </form>
            </p>
    </div>   
    </div>
</div>
</body>
</html>
