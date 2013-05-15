<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>SlimeBase How To</title>
<link href="../css/dbmain.css" rel="stylesheet" type="text/css" />
</head>

<body class="mainPage" onload="MM_preloadImages('images/layout/rollovers/home_over.jpg','images/layout/rollovers/background_over.jpg','images/layout/rollovers/pre_over.jpg','images/layout/rollovers/post_over.jpg','images/layout/rollovers/regulatory_over.jpg','images/layout/rollovers/proposal_over.jpg')">
<div id="container">
	<div id="dbHeader"></div>
    <?php include('../includes/dbmainNav.php');?>
	<div id="content" style="clear:both;">
    	<div id="memberLogIn"></div>
   	  <div id="leftContent">
   	    <p><a href='../database.php'><- return to main page</a><br />
   	      <br />
   	    </p>
      	<a name="top" id="top"></a>
        <h2 id="howToHeader">SlimeBase HowTo</h2>
        	<ol><li><a href="#browse">Browse Projects</a></li>
            <li><a href="#searchproj">Search Projects</a></li>
            <li><a href="#viewproj">View Project Sequences</a></li>
            <li><a href="#searchproj">Search Project Sequences</a></li>
            <li><a href="#viewseq">View a Sequence</a></li>
            <li><a href="#startcustom">Start a Custom FASTA</a></li>
            <li><a href="#editcustom">Edit a Custom FASTA</a></li>
            <li><a href="#customdownload">Download Custom FASTA</a></li>
            <li><a href="#customblast">Set Custom FASTA as BLAST Query</a></li>
            <li><a href="#performblast">Perform a BLAST</a></li>
            </ol>
            <br /><br />
   	    <h1><a name="browse" id="browse"></a>Browsing Projects</h1> 
		<p>Projects in SlimeBase are organized hierarchically into folders and are accessed via the project browser panel.        </p>
		<p>&nbsp;</p>
		<p><img src="tutorial_images/browse1.jpg" alt="" width="306" height="382" /><br />
	      <br />
	    To browse SlimeBase projects, click on the desired project name. This will open the project, displaying any 					additional projects contained within that project. Opening a project will also add it to the project path. </p>
		<p>&nbsp;</p>
		<p><img src="tutorial_images/browse2.jpg" alt="" width="285" height="309" /><br />
		  <br />
	    Clicking a project name on the project path will return you to that project. </p>
		<p><img src="tutorial_images/browse3.jpg" alt="" width="288" height="291" /><br />
	      <br />
	      <br />
	      <br />
	      <br />
	      <br />
	      <br />
	    </p>
		<h1><a name="searchproj" id="searchproj"></a>Searching a Project <a href="#top">return to top</a><br />
		  <br />
		</h1> 
        <p>To search a SlimeBase sequencing project, click the magnifying glass in the project browser.  </p>
        <p>&nbsp;</p>
        <p><img src="tutorial_images/projsearch1" alt="" /></p>
        <p>&nbsp;</p>
        <p>&nbsp;</p>
        <p>A search box will appear directly to the right of the magnifying glass.   Enter your search terms and hit the enter key to perform the search.<br />  
            <br />  
          <img src="tutorial_images/search2.jpg" alt="" width="288" height="291" /><br />
          <br />  
          </p>
        <p><br />
          Search results are listed below the project path. <br />
          <br />  
          <img src="tutorial_images/search3.jpg" alt="" width="288" height="291" /><br />
          <br />
          <br />
        To close the search result, click any project name in the project path or the home button. </p>
        <p>&nbsp;</p>
        <p><img src="tutorial_images/search4.jpg" alt="" width="288" height="291" /><br />
          <br />
          <br />
          <br />
          <br />
          <br />
        </p>
        <h1><a name="viewproj" id="viewproj"></a>Viewing a Project’s Sequences <a href="#top">return to top</a><br />
		</h1>
          <br />
		<p>To view project sequences, navigate to the project you wish to open in the project browser. If a project contains sequences, you will see a button named “Load Sequences” listed with the other contents of the project. Click this button to load the sequences into the sequence panel.<br />
	    </p><br />
		<p><img src="tutorial_images/projview1.jpg" alt="" width="288" height="291" /><br />
		  <br />
	    </p>
		<p>&nbsp;</p>
		<p>&nbsp;</p>
		<p>A sequence entry in SlimeBase has the following elements.</p>
		<p>&nbsp;</p>
	    <p><strong>a.</strong> Unique SlimeBase identifier</p>
	    <p><strong>b. </strong>Read/Contig name</p>
	    <p><strong>c.</strong> E-value from annotation attempt</p>
	    <p><strong>d. </strong>Annotation against NCBI's SP or NR database</p>
	    <p>&nbsp;</p>
	    <p><img src="tutorial_images/sequence_description.jpg" width="438" height="84" /><br />
		    <br />
	        <br />
		    <br />
	    </p>
      <h1><a name="searchseqs" id="searchseqs"></a>Searching a Project’s Sequences <a href="#top">return to top</a></h1>
		<p><br />
	    To search a project by keyword, first make sure the project is open in the sequence panel. Click the search link in the upper-right hand corner of the sequence panel to open the search dialogue box.  </p>
		<p>&nbsp;</p>
	  <p><br />
		    <img src="tutorial_images/seqsearch1.jpg" alt="" width="574" height="434" /><br />
	            <br />
              <br />
              Enter your search terms and either hit enter or the “Go” button.<br />
              <img src="tutorial_images/seqsearch2.jpg" alt="" width="574" height="434" /><br />
                  <br />
                  <br />
                Your search results will be listed in the sequence panel. To return to browsing the project, click the “return to project” link.<br />
                <br />
                <img src="tutorial_images/seqsearch3.jpg" alt="" width="574" height="434" /><br />
		          <br />
                  <br />
                  <br />
                  <br />
          <br />
	    </p>
		<h1><a name="viewseq" id="viewseq"></a>Viewing a Sequence <a href="#top">return to top</a></h1>
        <p><br />
          To view a sequence, you can either click the sequence identifier to load the sequence as its default type (NT or AA), or the NT or AA buttons to load the nucleotide or amino acid sequence, respectively. If SlimeBase does not contain an amino acid sequence for a given identifier, it will present the translation of the nucleotide sequence in all six reading frames. <br />
          <br />
          <img src="tutorial_images/seqview1.jpg" alt="" width="574" height="434" /><br />
          <br />
          <br />
          The sequence is loaded into the detail panel, which is directly below the sequence panel.<br />
          <br />
          <img src="tutorial_images/seqview2.jpg" alt="" width="561" height="584" /><br />
          <br />
          <br />
          <br />
          <br />
          <br />
        </p>
        <h1><a name="startcustom" id="startcustom"></a>Starting/Adding a Sequence to Custom FASTA <a href="#top">return to top</a><br />
	</h1>

<p><br />
  SlimeBase allows users to build and maintain a custom FASTA file as they browse sequences. This file can then be downloaded or set as the query to be used in a BLAST search. To start a custom FASTA, or add a sequence to an existing custom FASTA, click the (insert custom fasta button) next to the sequence identifier in the sequence panel. <br />
  <br />
  <img src="tutorial_images/startcustom1.jpg" alt="" width="566" height="437" /><br />
  <br />
  <br />
  The added sequence will appear in the Custom FASTA panel, which shares space with the sequence detail panel. To switch between custom FASTA and sequence detail views, click the link in the upper right hand side of the panel.<br />
  <br />
  <img src="tutorial_images/startcustom2.jpg" alt="" width="569" height="291" /><br />
  <br />
  <br />
  <br />
  <br />
</p>

<h1><a name="editcustom" id="editcustom"></a>Changing the Type/Removing a Sequence From a Custom FASTA <a href="#top">return to top</a></h1>

<p><br />
  To change the type of a custom FASTA sequence, select either NT or AA from the drop down menu in the custom FASTA window. This will change the type of the sequence for the purpose of downloading and performing a BLAST search.<br />
    <br />
        <img src="tutorial_images/editcustom1.jpg" alt="" width="569" height="291" /><br />
        <br /><br />

To remove a sequence a custom FASTA sequence, click the trash can icon next to the sequence identifier.<br />
<br />
<img src="tutorial_images/editcustom2.jpg" alt="" width="569" height="291" /><br />
<br />
<br />
<br />
<br />
<br />
</p>
<h1><a name="customdownload" id="customdownload"></a>Downloading a Custom FASTA <a href='#top'>return to top</a></h1>
<p><br />
  To download a custom FASTA file, click the download button in the Custom FASTA Panel. This will prepare  your sequences for download. Click the download link when the file is ready to download your sequences.<br />
  <br />
  <img src="tutorial_images/customdownload.jpg" alt="" width="566" height="250" /><br />
  <br />
  <br />
  <br />
</p>
<h1><a name="customblast" id="customblast"></a>Setting a Custom FASTA as a BLAST Query <a href='#top'>return to top</a><br />
  <br />
</h1>
<p>To set a custom FASTA as a query set for BLAST, click the set blast button in the Custom FASTA Panel. <br />
  <br />
  <img src="tutorial_images/setcustom.jpg" alt="" width="566" height="250" /><br />
  <br />
  The &quot;Use Custom FASTA&quot; button will check in the BLAST Panel. Remember that only one source of query sequences can be set at any given time.<br />
      <img src="tutorial_images/setcustom2.jpg" alt="" width="231" height="41" /><br />
          <br />
          <br />
          <br />
        <br />
  <br />
</p>
<h1><a name="performblast" id="performblast"></a>Performing a BLAST <a href="#top">return to top</a></h1>

<p><br />
  BLASTing at SlimeBase is performed via the blast panel, which is one the far right-hand side of the page. To perform an alignment, select the BLAST program, the desired database, the desired e-value, and set the desired source of query sequences. Query sequences at SlimeBase can either be pasted into the provided sequence box, uploaded in a FASTA file, or provided via a custom FASTA.  Up to 100 sequences can be manually pasted into the sequence box. The maximum size of an uploaded file is 10MB.
  <br />
  <br />
  <img src="tutorial_images/blast1.jpg" alt="" width="277" height="584" /><br />
  <br /><br />
If you provide an email address, SlimeBase will send you an email notifying you that your job has been completed, as well as provide a link to a zip file containing your results. Alignment results will automatically removed from our servers after three days.
<br /><br />
If no email is provided, you can wait at the SlimeBase page for your job to finish. Job status is displayed directly below the “Go” button on the blast panel, and gives your job’s current position in the queue. When your BLAST job is finished, you can click the status text to view your results within SlimeBase. 
<br />
<br />
<img src="tutorial_images/blast2.jpg" alt="" width="277" height="584" /><br />
<br /><br />
SlimeBases’s result visualization has three tabs, allowing you to view the textual BLAST output, the quantification data for your query sequences, and graphical representations of each alignment. Click the close button to return to SlimeBase’s main page.</p>
        <p>&nbsp;</p>
        <p>&nbsp;</p>
        <br />
   	  </div>
   	   <div class="quickLinks">
      	<div id="rightContent">
        <div id="linkContainer">
        <h1>Quick Links</h1>
      	  <p><a href="http://10.41.128.72/admin/seq_view/flash/tdb_prototype.html">SlimeBase 1.0</a><br />
          <a href="http://mascot.biotech.ufl.edu/bq/">BlastQuest</a>   	      </p>
      	  </div>
      	</div>
   	  </div>
  </div>
	<div id="footer"></div>
</div>
</body>
</html>
