<link rel="stylesheet" href="<?php echo $this->url () ?>/calendar.css" type="text/css" media="screen" title="no title" charset="utf-8"/>
<script type="text/javascript" src="<?php echo $this->url () ?>/js/microajax.js"></script>
<script type="text/javascript" charset="utf-8">
/*<![CDATA[ */
	function calendar (year,month,full,base)
	{
		microAjax ('<?php echo $this->url () ?>/ajax.php?full=' + full + '&month=' + month + '&year=' + year + '&base=' + base, function(response) { document.getElementById ('giraffe_calendar').innerHTML = response});
		return false;
	}
/*]]>*/
</script>
