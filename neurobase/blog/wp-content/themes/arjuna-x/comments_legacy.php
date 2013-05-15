<?php
// This is the comments file for Wordpress 2.6.x and older versions

// Forbid direct access
if (!empty($_SERVER['SCRIPT_FILENAME']) && 'comments_legacy.php' == basename($_SERVER['SCRIPT_FILENAME']))
	die ('This page cannot be loaded directly.');

// Password protection
if (!empty($post->post_password) && $_COOKIE['wp-postpass_' . COOKIEHASH] != $post->post_password) {
	?>
	<p class="noComments"><?php _e('This post is password protected. Enter the password to view comments.', 'Arjuna'); ?></p>
	<?php return;
}

//WP 2.6 and older Comment Loop
?>
<div class="commentHeader">
	<h4><?php _e('Comments', 'Arjuna'); ?></h4>
	<?php if(comments_open()): ?>
		<a href="#respond" class="btnReply"><span><?php _e('Leave a comment', 'Arjuna'); ?></span></a>
	<?php endif; ?>
	<?php if(pings_open()): ?>
		<a href="<?php trackback_url(); ?>" class="btnTrackback"><span><?php _e('Trackback', 'Arjuna'); ?></span></a>
	<?php endif; ?>
</div>
<?php if ($comments) { ?>
<ul class="commentList">
	<?php foreach ( $comments as $comment ): ?>
		<li class="comment <?php if($comment->comment_author_email == get_the_author_email()) {echo 'adminComment';} ?>" id="comment-<?php comment_ID() ?>">
			<div class="author">
				<div class="avatar">
					<?php 
					if (function_exists('get_avatar')) {
						echo get_avatar($comment, 60);
					} else {
						//gravatar code for < 2.5
						$gravUrl = "http://www.gravatar.com/avatar.php?gravatar_id=" . md5($email) . "&size=" . $size;
						echo "<img src='$gravUrl' height='60px' width='60px' />";
					 }
					?>
				</div>
				<div class="name">
					<?php if (get_comment_author_url()): ?>
						<a id="commentauthor-<?php comment_ID() ?>" class="url" href="<?php comment_author_url() ?>" rel="external nofollow">
					<?php else: ?>
						<span id="commentauthor-<?php comment_ID() ?>">
					<?php endif; ?>
					<?php comment_author(); ?>
					<?php if(get_comment_author_url()): ?>
						</a>
					<?php else: ?>
						</span>
					<?php endif; ?>
				</div>
			</div>
			<div class="messageBox">
				<div class="date">
					<?php printf( __('%1$s at %2$s', 'Arjuna'), get_comment_time(get_option('date_format')), get_comment_time(__('H:i', 'Arjuna')) ); ?>
				</div>
				<div class="links">
					<?php edit_comment_link('Edit','',''); ?>
				</div>
				<div class="content">
					<?php if ($comment->comment_approved == '0') : ?>
						<p><small><?php _e('Your comment is awaiting moderation.', 'Arjuna'); ?></small></p>
					<?php endif; ?>
	
					<div id="commentbody-<?php comment_ID() ?>">
						<?php comment_text(); ?>
					</div>
				</div>
			</div>
			
		</li>
	<?php endforeach; /* end for each comment */ ?>
</ul>

<?php // NOTE: NOT YET IMPLEMENTED PROPERLY INTO ARJUNA // ?>
<div class="commentNavigation">
	<?php if(function_exists('paginate_comments_links')) { ?>
		<?php paginate_comments_links('prev_text='.__('Previous', 'Arjuna').'&next_text='.__('Next', 'Arjuna').''); ?>
	<?php } else { ?>
		<div class="older"><?php previous_comments_link(__('Older Comments', 'Arjuna')) ?></div>
		<div class="newer"><?php next_comments_link(__('Newer Comments', 'Arjuna')) ?></div>
	<?php } ?>
</div>
<?php // NOTE END // ?>

<?php } else { // no comments (yet) ?>
	<?php if ('open' == $post->comment_status) { ?>
		<p class="noComments"><?php _e('No one has commented yet.', 'Arjuna'); ?></p>
	<?php } else { ?>
		<p class="noComments"><?php _e('Comments are closed.', 'Arjuna'); ?></p>
	<?php } ?>
<?php } ?>

<?php
// Commment form
if ('open' == $post->comment_status):
?>
<div class="commentReply" id="respond">
	<h4 class="replyHeader"><?php _e('Leave a Comment', 'Arjuna'); ?></h4>
	<?php if ( get_option('comment_registration') && !$user_ID ) : ?>
		<p style="margin-bottom:40px;"><?php printf(__('You must be %slogged in%s to post a comment.', 'Arjuna'), '<a href="'.get_option('siteurl').'/wp-login.php?redirect_to='.get_permalink().'">', '</a>'); ?></p></div>
	<?php else : ?>
		<form action="<?php echo get_option('siteurl'); ?>/wp-comments-post.php" name="reply" method="post" id="commentform">
			<input type="hidden" id="replyNameDefault" value="<?php _e('Your name', 'Arjuna'); ?>" />
			<input type="hidden" id="replyEmailDefault" value="<?php _e('Your email', 'Arjuna'); ?>" />
			<input type="hidden" id="replyURLDefault" value="<?php _e('Your website', 'Arjuna'); ?>" />
			<input type="hidden" id="replyMsgDefault" value="<?php _e('Your comment', 'Arjuna'); ?>" />
			<?php if ( $user_ID ): ?>
			<p><?php printf(__('Logged in as %s.', 'Arjuna'), '<a href="'.get_option('siteurl').'/wp-admin/profile.php">'.$user_identity.'</a>'); ?> <a href="<?php echo wp_logout_url(get_permalink()); ?>" title="<?php _e('Log out of this account', 'Arjuna'); ?>"><?php _e('Logout', 'Arjuna'); ?> &raquo;</a></p>
			<?php else : ?>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author)): ?> inputIA<?php endif; ?>" id="replyName" name="author" value="<?php if(!empty($comment_author)) { echo $comment_author; } else { _e('Your name', 'Arjuna'); } ?>" /></div>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author_email)): ?> inputIA<?php endif; ?>" id="replyEmail" name="email" value="<?php if(!empty($comment_author)) { echo $comment_author_email; } else { _e('Your email', 'Arjuna'); } ?>" /></div>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author_url)): ?> inputIA<?php endif; ?>" id="replyURL" name="url" value="<?php if(!empty($comment_author_url)) { echo $comment_author_url; } else { _e('Your website', 'Arjuna'); } ?>" /></div>
			<?php endif; ?>
			<input type="hidden" name="comment_post_ID" value="<?php echo $id; ?>" />
			<div class="replyRow"><textarea class="inputIA" id="replyMsg" name="comment"><?php _e('Your comment', 'Arjuna'); ?></textarea></div>
			<div class="replySubmitArea">
				<a href="<?php echo get_post_comments_feed_link(); ?>" class="btnSubscribe"><span><?php _e('Subscribe to comments', 'Arjuna'); ?></span></a>
				<button type="submit" class="inputBtn" value="Submit" name="submit"><?php _e('Leave comment', 'Arjuna'); ?></button>
			</div>
			<?php do_action('comment_form', $post->ID); ?>
		</form>
	<?php endif; // If registration required and not logged in ?>
</div>
<?php endif; ?>
