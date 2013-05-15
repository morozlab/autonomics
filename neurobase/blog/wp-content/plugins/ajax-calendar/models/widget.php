<?php

class AJAX_Calendar_Widget extends Widget_AJAX_Calendar
{
	var $title = '';
	var $category_ids = '';

	function has_config () { return true; }

	function load ($config)
	{
		if (isset ($config['title']))
			$this->title = $config['title'];
		if (isset ($config['category_id']))
			$this->category_ids = $config['category_id'];
	}

	function display ($args)
	{
		extract ($args);

		echo $before_widget;

		if ($this->title)
			echo $before_title.$this->title.$after_title;

		ajax_calendar ($this->category_ids);

		echo $after_widget;
	}

	function description ()
	{
		return 'AJAX-powered calendar';
	}

	function config ($config, $pos)
	{
		?>
		<table>
			<tr>
				<th>Title:</th>
				<td><input type="text" name="<?php echo $this->config_name ('title', $pos) ?>" value="<?php echo htmlspecialchars ($config['title']) ?>"/></td>
			</tr>
			<tr>
				<th>Category ID:</th>
				<td><input type="text" name="<?php echo $this->config_name ('category_id', $pos) ?>" value="<?php echo htmlspecialchars ($config['category_id']) ?>"/></td>
			</tr>
		</table>
		<?php
	}

	function save ($data)
	{
		return array ('title' => $data['title'], 'category_id' => $data['category_id']);
	}
}

?>