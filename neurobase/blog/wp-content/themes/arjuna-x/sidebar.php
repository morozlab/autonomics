<?php $arjunaOptions = arjuna_get_options(); ?>
<?php
if ($arjunaOptions['sidebarDisplay']!='none') {
?><div class="sidebars">
	<div class="t"><div></div></div>
	<div class="i"><div class="i2"><div class="c">
		<div class="sidebarIcons">
			<?php if($arjunaOptions['sidebar_showRSSButton']): ?><a class="rssBtn" href="<?php bloginfo('rss2_url'); ?>"><?php if($arjunaOptions['sidebar_displayButtonTexts']): ?>RSS<?php endif; ?></a><?php endif; ?>
			<?php if($arjunaOptions['sidebar_showTwitterButton']): ?>
			<a class="twitterBtn" href="<?php echo $arjunaOptions['sidebar_twitterURL']; ?>"><?php if($arjunaOptions['sidebar_displayButtonTexts']): ?>Twitter<?php endif; ?></a>
			<?php endif; ?>
			<?php if($arjunaOptions['sidebar_showFacebookButton']): ?>
			<a class="facebookBtn" href="<?php echo $arjunaOptions['sidebar_facebookURL']; ?>"><?php if($arjunaOptions['sidebar_displayButtonTexts']): ?>Facebook<?php endif; ?></a>
			<?php endif; ?>
		</div>
		<div>
		<?php if (!function_exists('dynamic_sidebar') || !dynamic_sidebar('Sidebar Top')): ?>
			<?php if($arjunaOptions['sidebar_showDefault']): ?>
			<div class="sidebarbox">
			<h4><span><?php _e('Recent Posts', 'Arjuna'); ?></span></h4>
			<ul>
			<?php wp_get_archives('type=postbypost&limit=10'); ?>
			</ul>
			</div>
			
			<div class="sidebarbox">
			<h4><span><?php _e('Browse by Tags', 'Arjuna'); ?></span></h4>
			<?php wp_tag_cloud('smallest=8&largest=17&number=30'); ?>
			</div>
			<?php endif; ?>
		<?php endif; ?>
		</div>
		<div class="sidebarLeft">
		<?php if (!function_exists('dynamic_sidebar') || !dynamic_sidebar('Sidebar Left')): ?>
			<?php if($arjunaOptions['sidebar_showDefault']): ?>
			<div class="sidebarbox">
			<h4><span><?php _e('Categories', 'Arjuna'); ?></span></h4>
			<ul>
				<?php wp_list_categories('show_count=0&title_li='); ?>
			</ul>
			</div>
			<?php endif; ?>
		<?php endif; ?>
		</div>
		<div class="sidebarRight">
		<?php if (!function_exists('dynamic_sidebar') || !dynamic_sidebar('Sidebar Right')): ?>
			<?php if($arjunaOptions['sidebar_showDefault']): ?>
			<div class="sidebarbox">
			<h4><span><?php _e('Meta', 'Arjuna'); ?></span></h4>
			<ul>
				<?php wp_register(); ?>
				<li><?php wp_loginout(); ?></li>
				<li><a href="http://validator.w3.org/check/referer" title="This page validates as XHTML 1.0 Transitional">Valid <abbr title="eXtensible HyperText Markup Language">XHTML</abbr></a></li>
				<?php wp_meta(); ?>
			</ul>
			</div>
			<?php endif; ?>
		<?php endif; ?>
		</div>
		<div class="clear">
			<?php if ( !function_exists('dynamic_sidebar') || !dynamic_sidebar('Sidebar Bottom') ) : ?>
			<?php endif; ?>
		</div>
	</div></div></div>
	<div class="b"><div></div></div>
</div><?php
}
?>

