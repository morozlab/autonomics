<?php $arjunaOptions = arjuna_get_options(); ?>
		<div class="clear"></div>
	</div><!-- .contentWrapper -->
	<div class="<?php if($arjunaOptions['footerStyle']=='style1'): ?>footer<?php else: ?>footer2<?php endif; ?>">
		<a href="http://www.wordpress.org" class="icon1"><img src="<?php bloginfo('template_url'); ?>/images/<?php if($arjunaOptions['footerStyle']=='style1'): ?>wordpressIcon.png<?php else: ?>wordpressIcon2.jpg<?php endif; ?>" width="20" height="20" alt="Powered by WordPress" /></a>
		<a class="icon2"><img src="<?php bloginfo('template_url'); ?>/images/<?php if($arjunaOptions['footerStyle']=='style1'): ?>srsIcon.png<?php else: ?>srsIcon2.jpg<?php endif; ?>" width="18" height="30" alt="Web Design by SRS Solutions" /></a>
		<span class="copyright">&copy; <?php print date('Y'); ?> <?php bloginfo('name'); ?></span>
		<span class="design"><a href="http://www.srssolutions.com/en/" title="Web Design by SRS Solutions">Design by <em>SRS Solutions</em></a></span>
	</div>
</div><!-- .pageContainer -->

<?php wp_footer(); ?>
</body>
</html>
