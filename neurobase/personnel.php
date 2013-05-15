<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Lab Personnel</title>
<link href="css/main.css" rel="stylesheet" type="text/css" />
<script type="text/javascript">
<!--
function MM_swapImgRestore() { //v3.0
  var i,x,a=document.MM_sr; for(i=0;a&&i<a.length&&(x=a[i])&&x.oSrc;i++) x.src=x.oSrc;
}
function MM_preloadImages() { //v3.0
  var d=document; if(d.images){ if(!d.MM_p) d.MM_p=new Array();
    var i,j=d.MM_p.length,a=MM_preloadImages.arguments; for(i=0; i<a.length; i++)
    if (a[i].indexOf("#")!=0){ d.MM_p[j]=new Image; d.MM_p[j++].src=a[i];}}
}

function MM_findObj(n, d) { //v4.01
  var p,i,x;  if(!d) d=document; if((p=n.indexOf("?"))>0&&parent.frames.length) {
    d=parent.frames[n.substring(p+1)].document; n=n.substring(0,p);}
  if(!(x=d[n])&&d.all) x=d.all[n]; for (i=0;!x&&i<d.forms.length;i++) x=d.forms[i][n];
  for(i=0;!x&&d.layers&&i<d.layers.length;i++) x=MM_findObj(n,d.layers[i].document);
  if(!x && d.getElementById) x=d.getElementById(n); return x;
}

function MM_swapImage() { //v3.0
  var i,j=0,x,a=MM_swapImage.arguments; document.MM_sr=new Array; for(i=0;i<(a.length-2);i+=3)
   if ((x=MM_findObj(a[i]))!=null){document.MM_sr[j++]=x; if(!x.oSrc) x.oSrc=x.src; x.src=a[i+2];}
}
//-->
</script>
</head>

<body class="mainPage" onload="MM_preloadImages('images/layout/rollovers/home_over.jpg','images/layout/rollovers/background_over.jpg','images/layout/rollovers/pre_over.jpg','images/layout/rollovers/post_over.jpg','images/layout/rollovers/regulatory_over.jpg','images/layout/rollovers/proposal_over.jpg','images/layout/rollovers/moroz_lab_member_roll.jpg')">
<div id="container">
	<div id="header"></div>
    <?php include('includes/mainNav.php');?>
	<div id="content" style="clear:both;">
    	<div id="memberLogIn"><a href="lab_members/member_log_in.html" onmouseout="MM_swapImgRestore()" onmouseover="MM_swapImage('Member_Roll','','images/layout/rollovers/moroz_lab_member_roll.jpg',1)"><img src="images/layout/rollovers/moroz_lab_member.jpg" alt="Click to log in" name="Member_Roll" width="144" height="18" border="0" id="Member_Roll" /></a></div>
   	  <div id="leftContent">
   	    <h1>Personnel</h1>
        <p><a href="mailto:moroz@whitney.ufl.edu">Leonid L. Moroz</a>, Ph.D. – Professor of Neuroscience<br />
		  <a href="mailto:jcn@whitney.ufl.edu">James C. Netherton III</a>, Chemist – chemistry, lab management<br />
		  </p>
        <p>Yelena Bobkova – Biological Scientist</p>
        <p><a href="mailto:fyu@ufl.edu">Fahong Yu</a>, Bioinformatics and Software Specialist – bioinformatics, database<br />
		Do Sung Sohn, Graduate Student<br />
        </p>
        <h3>&nbsp;</h3>
        <h3>
	    Matt Citarella, Graduate Student – computer information science and engineering</h3><p><br />
	      <img id="lFloatImg" src="images/personnel/matt_bot.jpg" alt="Matt Citarella" width="175" height="130" />	      Matt is a second year masters student in the CISE Bioinformatics program
at the University of Florida. He received his Bachelor of Arts in
Digital Media from the University of Central Florida in 2006, the same
year he joined the Moroz Lab. Matt performs the lions share of the lab's
Bioinformatics work including: phylogenetic analysis, digital gene
expression profiling, sequence annotation, sequence alignment, database
curation, web design, and web application development.

            <br />
            <br />
          Mat's primary areas of interest are human-usable interface design, data
presentation techniques, gene analysis, web development, food, and
sleep. He is currently designing a web-based transcriptome sequence and
annotation viewer to handle the products of the lab's various sequencing
projects.<br />
<br />
email: <a href="mailto:mcitarel@cise.ufl.edu">mcitarel@cise.ufl.edu</a><br />
</p>
	    <p>phone: 904-461-4007
          <br />
		</p>
	    <p>&nbsp;</p>
	    </h3>
	  </div>
   	  <div class="quickLinks">
   	    <div id="rightContent"></div>
   	  </div>
  </div>
	<div id="footer"></div>
</div>
</body>
</html>
