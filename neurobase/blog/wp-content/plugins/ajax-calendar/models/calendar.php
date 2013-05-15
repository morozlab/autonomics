<?php

class Calendar
{
	var $categories = array();
	var $url;
	var $showit = true;
	var $split_open  = "&laquo;-&raquo;";
	var $split_close = "&raquo;-&laquo;";

	function Calendar ($url, $categories = '' )
	{
		$this->url = $url;

		if ($categories)
			$this->categories = explode (',', $categories);
	}

	function show ($month_num = '', $year_num = '', $full = '0', $categories = '')
  {
    global $wpdb, $m, $monthnum, $year, $timedifference, $month, $month_abbrev, $weekday, $weekday_initial, $weekday_abbrev, $posts, $wp_rewrite;

    $month_num = intval ($month_num);
    $year_num  = intval ($year_num);

    if ($year_num <= 2003)
      $year_num = date ('Y');

    if ($month_num <= 0 || $month_num > 12)
      $month_num = date ('n');

    $monthnum = $month_num;
    $year     = $year_num;

    // Is this cached?
		$cache = wp_cache_get ($year."_".$monthnum."_".$full, 'calendar');
    if ($cache)
			$text = $cache;
    else
    {
      $next_month  = $monthnum + 1;
      $next_year   = $year;
      $prev_month  = $monthnum - 1;
      $prev_year   = $year;
      $notfull     = $full == 0 ? 1 : 0;

      if ($next_month == 13)
      {
        $next_month = 1;
        $next_year++;
      }

      if ($prev_month == 0)
      {
        $prev_month = 12;
        $prev_year--;
      }

      ob_start ();
      $this->wp_calendar ($month_num, $year_num, $full);
      $text = ob_get_contents ();
      ob_end_clean ();

      // Cache this
			wp_cache_add ($year."_".$monthnum."_".$full, $text, 'calendar');
    }

    return $text;
  }


  // This is a duplicate of the WP function, but modified slightly to provide the enhanced functionality
  function wp_calendar ($thismonth, $thisyear, $full)
  {
  	global $wpdb, $m, $monthnum, $year, $timedifference, $month, $month_abbrev, $weekday, $weekday_initial, $weekday_abbrev, $posts;
		global $giraffe_ajax_always;

		if ($giraffe_ajax_always && $_SERVER['HTTP_X_REQUESTED_WITH'] != 'XMLHttpRequest')
			$full = 1;

  	// week_begins = 0 stands for Sunday
  	$week_begins = intval(get_option('start_of_week'));
  	$add_hours = intval(get_option('gmt_offset'));
  	$add_minutes = intval(60 * (get_option('gmt_offset') - $add_hours));

  	// Let's figure out when we are
  	$unixmonth = mktime(0, 0 , 0, $thismonth, 1, $thisyear);

	// Get the current post
    $current = $wpdb->get_row ("SELECT count(post_date) as count FROM $wpdb->posts WHERE post_date='".date('Y')."-".date ('n')."-01' AND post_status='publish'  AND post_type='post'");

  	// Get the next and previous month and year with at least one post
  	$previous = $wpdb->get_row("SELECT DISTINCT MONTH(post_date) AS month, YEAR(post_date) AS year
  		FROM $wpdb->posts
  		WHERE post_date < '$thisyear-$thismonth-01'
  		AND post_status = 'publish' AND post_type='post'
  			ORDER BY post_date DESC
  			LIMIT 1");
  	$next = $wpdb->get_row("SELECT	DISTINCT MONTH(post_date) AS month, YEAR(post_date) AS year
  		FROM $wpdb->posts
  		WHERE post_date >	'$thisyear-$thismonth-01'
  		AND MONTH( post_date ) != MONTH( '$thisyear-$thismonth-01' )
  		AND post_status = 'publish' AND post_type='post'
  			ORDER	BY post_date ASC
  			LIMIT 1");

    if ($current->count == 0)
      $initiallink = $month[zeroise($thismonth, 2)].' '.$thisyear;
    else
      $initiallink = '<a href="'.get_month_link (date ('Y'), date ('n')).'" onclick="calendar ('.date ('n').','.date ('Y').',\''.$full.'\');return false;">' . $month[zeroise($thismonth, 2)] . ' ' . date('Y', $unixmonth) . '</a>';

  	echo '<table id="wp-calendar">
  	<caption id="wp-calendar-caption">'.$initiallink.'</caption>
  	<thead>
  	<tr>';

  	$day_abbrev = $weekday_initial;
  	// if ( $daylength > 1 )
  	// 	$day_abbrev = $weekday_abbrev;

  	$myweek = array();

  	for ( $wdcount=0; $wdcount<=6; $wdcount++ ) {
  		$myweek[]=$weekday[($wdcount+$week_begins)%7];
  	}

  	foreach ( $myweek as $wd ) {
  		echo "\n\t\t<th abbr=\"$wd\" scope=\"col\" title=\"$wd\">" . $day_abbrev[$wd] . '</th>';
  	}

  	echo '</tr></thead><tfoot><tr>';

  	if ($previous)
  	{
  	  $link = get_month_link($previous->year, $previous->month);
  	  $linktext = sprintf(__('View posts for %1$s %2$s'), $month[zeroise($previous->month, 2)], date('Y', mktime(0, 0 , 0, $previous->month, 1, $previous->year)));
   	  ?>
   		<td abbr="<?php $month[zeroise($previous->month, 2)]?>" colspan="3" id="prev">
   		  <a href="<?php echo $link ?>" title="<?php echo $linktext ?>" onclick="return calendar(<?php echo $previous->year ?>,<?php echo $previous->month ?>,<?php echo $full ?>,'<?php echo implode( ',', $this->categories ); ?>')">
   		  &laquo; <?php echo $month_abbrev[$month[zeroise($previous->month, 2)]] ?>
   		  </a></td>
   		<?php
   	}
   	else

  		echo "\n\t\t".'<td colspan="3" id="prev" class="pad">&nbsp;</td>';

	if ($this->showit)
	{
		if ( !is_array( $this->categories ) )
			$this->categories = array();
     echo '<td id="showit"><a href="'.get_month_link ($thisyear, $thismonth).'" onclick="return calendar ('.$thisyear.','.$thismonth.','.($full == 0 ? 1 : 0).',\''.implode( ',', $this->categories ).'\')">'.$this->split_open.'</a></td>';

  	if ( $next ) {
  		echo "\n\t\t".'<td abbr="' . $month[zeroise($next->month, 2)] . '" colspan="3" id="next"><a href="' .
  		get_month_link($next->year, $next->month) . '" title="View posts for ' . $month[zeroise($next->month, 2)] . ' ' .
  		date('Y', mktime(0, 0 , 0, $next->month, 1, $next->year)) .
  		'" onclick="return calendar ('.$next->year.','.$next->month.','.$full.')"'.
  		'>' . $month_abbrev[$month[zeroise($next->month, 2)]] . ' &raquo;</a></td>';
  	} else {
  		echo "\n\t\t".'<td colspan="3" id="next" class="pad">&nbsp;</td>';
  	}
	}
	
  	echo '
  	</tr>
  	</tfoot>

  	<tbody>
  	<tr>';

  	// Get days with posts
	$sql = "SELECT DISTINCT DAYOFMONTH(post_date)
  		FROM $wpdb->posts 
		LEFT JOIN {$wpdb->prefix}term_relationships ON ({$wpdb->prefix}posts.ID = {$wpdb->prefix}term_relationships.object_id)
		LEFT JOIN {$wpdb->prefix}term_taxonomy ON ({$wpdb->prefix}term_relationships.term_taxonomy_id = {$wpdb->prefix}term_taxonomy.term_taxonomy_id)
  		WHERE MONTH(post_date) = $thismonth
  		AND YEAR(post_date) = $thisyear
  		AND post_status = 'publish' AND post_type='post'";
		
	if (count ($this->categories) > 0)
		$sql .= " AND {$wpdb->prefix}term_taxonomy.term_id IN (".implode (',', $this->categories).')';

  	$dayswithposts = $wpdb->get_results($sql, ARRAY_N);
		
  	if ( $dayswithposts ) {
  		foreach ( $dayswithposts as $daywith ) {
  			$daywithpost[] = $daywith[0];
  		}
  	} else {
  		$daywithpost = array();
  	}



  	if ( strstr($_SERVER['HTTP_USER_AGENT'], 'MSIE') || strstr(strtolower($_SERVER['HTTP_USER_AGENT']), 'camino') || strstr(strtolower($_SERVER['HTTP_USER_AGENT']), 'safari') )
  		$ak_title_separator = "\n";
  	else
  		$ak_title_separator = ', ';

  	$ak_titles_for_day = array();
  	$ak_post_titles = $wpdb->get_results("SELECT post_title, DAYOFMONTH(post_date) as dom "
  		."FROM $wpdb->posts "
  		."WHERE YEAR(post_date) = '$thisyear' "
  		."AND MONTH(post_date) = '$thismonth' "
  		."AND post_status = 'publish' AND post_type='post'"
  	);
  	if ( $ak_post_titles ) {
  		foreach ( $ak_post_titles as $ak_post_title ) {
  				if ( empty($ak_titles_for_day['day_'.$ak_post_title->dom]) )
  					$ak_titles_for_day['day_'.$ak_post_title->dom] = '';
  				if ( empty($ak_titles_for_day["$ak_post_title->dom"]) ) // first one
  					$ak_titles_for_day["$ak_post_title->dom"] = str_replace('"', '&quot;', wptexturize($ak_post_title->post_title));
  				else
  					$ak_titles_for_day["$ak_post_title->dom"] .= $ak_title_separator . str_replace('"', '&quot;', wptexturize($ak_post_title->post_title));
  		}
  	}


  	// See how much we should pad in the beginning
  	$pad = calendar_week_mod(date('w', $unixmonth)-$week_begins);
  	if ( 0 != $pad )
  		echo "\n\t\t".'<td colspan="'.$pad.'" class="pad">&nbsp;</td>';

  	$daysinmonth = intval(date('t', $unixmonth));
  	for ( $day = 1; $day <= $daysinmonth; ++$day ) {
  		if ( isset($newrow) && $newrow )
  			echo "\n\t</tr>\n\t<tr>\n\t\t";
  		$newrow = false;

  		if ( $day == gmdate('j', (time() + (get_option('gmt_offset') * 3600))) && $thismonth == gmdate('m', time()+(get_option('gmt_offset') * 3600)) && $thisyear == gmdate('Y', time()+(get_option('gmt_offset') * 3600)) )
  			echo '<td id="today">';
  		else
  			echo '<td>';

  		if ( in_array($day, $daywithpost) ) // any posts today?
  				echo '<a href="' . get_day_link($thisyear, $thismonth, $day) . "\" title=\"$ak_titles_for_day[$day]\">$day</a>";
  		else
  			echo $day;
  		echo '</td>';

  		if ( 6 == calendar_week_mod(date('w', mktime(0, 0 , 0, $thismonth, $day, $thisyear))-$week_begins) )
  			$newrow = true;
  	}

  	$pad = 7 - calendar_week_mod(date('w', mktime(0, 0 , 0, $thismonth, $day, $thisyear))-$week_begins);
  	if ( $pad != 0 && $pad != 7 )
  		echo "\n\t\t".'<td class="pad" colspan="'.$pad.'">&nbsp;</td>';

  	echo "\n\t</tr>";
  	echo "\n\t</tbody>\n\t</table>";

    if ($full == 1)
    {
      global $wpdb;

      $res = $wpdb->get_results ("SELECT ID,post_title FROM $wpdb->posts WHERE month(post_date) = '$monthnum' AND year(post_date) = '$year' AND post_status = 'publish'  AND post_type='post' ORDER BY post_date");
      if ($res)
      {
        $text .= "<div id=\"wp-calendar-split\"><ul>";
        foreach ($res AS $pres)
          $text .= "<li><a href=\"".get_permalink ($pres->ID)."\">".apply_filters ('the_title', $pres->post_title)."</a></li>";

        $text .= "</ul></div>";
      }

      echo $text;
    }
  }
}
?>