<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head><?php
include("../includes/checkSession.php");
?>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Administrator Home</title>

<script type="text/javascript">
<!--
function MM_preloadImages() { //v3.0
  var d=document; if(d.images){ if(!d.MM_p) d.MM_p=new Array();
    var i,j=d.MM_p.length,a=MM_preloadImages.arguments; for(i=0; i<a.length; i++)
    if (a[i].indexOf("#")!=0){ d.MM_p[j]=new Image; d.MM_p[j++].src=a[i];}}
}
//-->
</script>
<link href="../css/main.css" rel="stylesheet" type="text/css" />
</head>

<body class="mainPage" onload="MM_preloadImages('../images/layout/rollovers/home_over.jpg','../images/layout/rollovers/background_over.jpg','../images/layout/rollovers/pre_over.jpg','../images/layout/rollovers/post_over.jpg','../images/layout/rollovers/regulatory_over.jpg','../images/layout/rollovers/proposal_over.jpg','../images/layout/rollovers/moroz_lab_member_roll.jpg')">
<div id="container">
  <div id="header"></div>
   <?php include('includes/admin_nav.php');?>
<div id="content">
      <div id="leftContent"> 
        <h1>Things to add...</h1>
        <p>Usage statistics<br />
        System Status<br />
        Hard Drive Space Left<br />
        Size of Database<br />
        Number of Sequences in Database</p>
      </div>
   	  <div id="rightContent">
   	</div>
  </div>
	<div id="footer"></div>
</div>
</body>
</html>