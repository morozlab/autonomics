<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>SlimeBase</title>
<script language='JavaScript'>
function getFlashVersion(){ 
  // ie 
  try { 
    try { 
      // avoid fp6 minor version lookup issues 
      // see: http://blog.deconcept.com/2006/01/11/getvariable-setvariable-crash-internet-explorer-flash-6/ 
      var axo = new ActiveXObject('ShockwaveFlash.ShockwaveFlash.6'); 
      try { axo.AllowScriptAccess = 'always'; } 
      catch(e) { return '6,0,0'; } 
    } catch(e) {} 
    return new ActiveXObject('ShockwaveFlash.ShockwaveFlash').GetVariable('$version').replace(/\D+/g, ',').match(/^,?(.+),?$/)[1]; 
  // other browsers 
  } catch(e) { 
    try { 
      if(navigator.mimeTypes["application/x-shockwave-flash"].enabledPlugin){ 
        return (navigator.plugins["Shockwave Flash 2.0"] || navigator.plugins["Shockwave Flash"]).description.replace(/\D+/g, ",").match(/^,?(.+),?$/)[1]; 
      } 
    } catch(e) {} 
  } 
  return '0,0,0'; 
} 
 
var version = getFlashVersion().split(',').shift(); 
if(version < 10){ 
  window.location = "getFlash.html";
}
</script>
<script src="../Scripts/AC_RunActiveContent.js" type="text/javascript"></script>
</head>

<body style="padding: 0px; margin:0px;">
<div align="center"><?php $ip=$_SERVER['REMOTE_ADDR'];
?>

  <script type="text/javascript">
AC_FL_RunContent( 'codebase','http://download.macromedia.com/pub/shockwave/cabs/flash/swflash.cab#version=9,0,28,0','width','1180','height','1200','src','../admin/seq_view/flash/tdb_prototype','flashvars','ip=<?php echo($ip); ?>','quality','high','pluginspage','http://www.adobe.com/shockwave/download/download.cgi?P1_Prod_Version=ShockwaveFlash','movie','../admin/seq_view/flash/tdb_prototype' ); //end AC code
</script>
  <noscript>
  Viewing SlimeBase requires that you have JavaScript enabled. Please enable it in your browser and try again.
  </noscript>
</div>
</body>
</html>
