<?php

include (dirname (__FILE__).'/../../../wp-config.php');

if (!class_exists ('AJAX_Calendar'))
	include (ABSPATH.'wp-content/plugins/ajax-calendar/ajax-calendar.php');
	
include (dirname (__FILE__).'/models/calendar.php');

$calendar_plugin = AJAX_Calendar::get ();
$cats = isset( $_GET['base'] ) ? $_GET['base'] : '';

$calendar = new Calendar ($calendar_plugin->url (), $cats );
echo $calendar->show (intval ($_GET['month']), intval ($_GET['year']), intval ($_GET['full']));

?>