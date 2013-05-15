<?php $arjunaOptions = arjuna_get_options(); ?>
<?php get_header(); ?>

<div class="contentArea">
	<?php if (have_posts()) : ?>
	<?php while (have_posts()) : the_post(); ?>
	<?php if($arjunaOptions['posts_showTopPostLinks'] && arjuna_has_other_posts()): ?>
	<div class="pagination paginationTop"><div>
		<?php arjuna_get_previous_post_link(__('Previous Post', 'Arjuna')); ?>
		<?php arjuna_get_next_post_link(__('Next Post', 'Arjuna')); ?>
	</div></div>
	<?php endif; ?>
	
	<div class="post" id="post-<?php the_ID(); ?>">
		<div class="postHeader">
			<h1 class="postTitle"><span><a href="<?php the_permalink() ?>" title="<?php _e('Permalink to', 'Arjuna'); ?> <?php the_title(); ?>"><?php the_title(); ?></a></span></h1>
			<div class="bottom"><div>
				<span class="postDate"><?php the_time(get_option('date_format')); ?><?php
					//Time
					if($arjunaOptions['postsShowTime']) {
						print _e(' at ', 'Arjuna'); the_time(get_option('time_format'));
					}
				?></span>
				<?php if($arjunaOptions['postsShowAuthor']): ?>
				<span class="postAuthor"><?php the_author_posts_link(); ?></span>
				<?php endif; ?>
				<?php if(!$arjunaOptions['comments_hideWhenDisabledOnPosts'] || ( 0 != $post->comment_count || comments_open() || pings_open() )): ?>
				<a href="<?php comments_link(); ?>" class="postCommentLabel"><span><?php
					if (function_exists('post_password_required') && post_password_required()) {
						_e('Pass required', 'Arjuna');
					} elseif(0 == $post->comment_count && !comments_open() && !pings_open()) {
						_e('Comments off', 'Arjuna'); 
					} else {
						comments_number(__('No comments', 'Arjuna'), __('1 comment', 'Arjuna'), __('% comments', 'Arjuna'));
					}
				?></span></a>
				<?php endif; ?>
			</div></div>
		</div>
		<div class="postContent">
			<?php the_content(__('continue reading...', 'Arjuna')); ?>
		</div>
		<div class="postFooter"><div class="r"></div>
			<div class="left">
				<span class="postCategories"><?php the_category(', '); ?></span>
				<?php if ( function_exists('the_tags') ): ?>
				<span class="postTags"><?php
					if (get_the_tags()) the_tags('', ', ', '');
					else print '<span>'.__('<i>none</i>', 'Arjuna').'</span>';
				?></span>
				<?php endif; ?>
			</div>
			<?php print arjuna_get_edit_link(__('Edit in Admin', 'Arjuna')); ?>
		</div>
	</div>
	<?php if(!$arjunaOptions['comments_hideWhenDisabledOnPosts'] || ( 0 != $post->comment_count || comments_open() || pings_open() )): ?>
	<div class="postComments" id="comments">
		<?php comments_template(); ?>
	</div>
	<?php endif; ?>

	<?php if($arjunaOptions['posts_showBottomPostLinks'] && arjuna_has_other_posts()): ?>
	<div class="pagination"><div>
		<?php arjuna_get_previous_post_link(__('Previous Post', 'Arjuna')); ?>
		<?php arjuna_get_next_post_link(__('Next Post', 'Arjuna')); ?>
	</div></div>
	<?php endif; ?>
	<?php endwhile; ?>



	<?php else : ?>
  <p><?php _e('There is nothing here.', 'Arjuna'); ?></p>
	<?php endif; ?>
</div>

<?php get_sidebar(); ?>

<?php get_footer(); ?>
