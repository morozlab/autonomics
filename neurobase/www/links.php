<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>Moroz Lab Links</title>
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
   	  <div id="leftContent" style="width:802px">
      	<p style="margin:0px;">Below you will find some links frequently used by our lab. Enjoy.</p>
        <h1>Aplysia Genome </h1>
        	<p>
            <a href="http://www.ncbi.nlm.nih.gov/sites/entrez?Db=genomeprj&cmd=ShowDetailView&TermToSearch=13634" target="_blank">Aplysia Genome at NCBI</a><br />
        	<a href="http://www.broad.mit.edu/node/435" target="_blank">Aplysia Genome at the Broad Institute</a><br />
            <a href="http://aplysia.cu-genome.org/html/index.html" target="_parent">Aplysia Genome Collaboration Page</a><br />
            <a href="http://www.genome.gov/13014443" target="_blank">NIH Approval of Aplysia Genome Project</a><br />
            <a href="http://www.genome.gov/10002154" target="_blank">NIH Approved Sequencing Targets</a><br />
		</p><br />
   	    <h1>Aplysia EST </h1>
        	<p>
            <a href="publications.php">Publications</a><br />
            <a href="http://mascot.biotech.ufl.edu/bq/" target="_blank">Aplysia Transcriptome @ UF's Genomic Core</a><br />
            <br />
        <h1>Aplysia Mitochondrial Genome</h1>
        	<p>
            <a href="publications.php">Publications</a><br />
            </p>
        <h1><br />
        Other Links</h1>
   	  <p>
            <a href="http://www.seaslugforum.net/factsheet.cfm?base=aplycali" target='_blank'>Slug Forum</a> <br />        
   	  		<a href="http://aplysia.miami.edu/" target="_blank">The National Resource for Aplysia</a><br />
	    <a href="http://www.jgi.doe.gov/">Joint Genome Institute</a><br />
        <a href="http://www.med.ufl.edu/IDP/">IDP UF</a><br />
<a href="http://www.ncbi.nlm.nih.gov/pubmed/">Pubmed</a><br />
<a href="http://www.ncbi.nlm.nih.gov/">NCBI</a><br />
 	    </p>
   	  </div>
   	  </div>
	<div id="footer"></div>
</div>
</body>
</html>
