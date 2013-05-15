<?php $arjunaOptions = arjuna_get_options(); ?>
<?php
// This is the comments file for Wordpress 2.7+

// Forbid direct access
if (!empty($_SERVER['SCRIPT_FILENAME']) && 'comments.php' == basename($_SERVER['SCRIPT_FILENAME']))
	die ('This page cannot be loaded directly.');

// Password protection
if (function_exists('post_password_required')) {
	if ( post_password_required() ) {
		?>
		<p class="noComments"><?php _e('This post is password protected. Enter the password to view comments.', 'Arjuna'); ?></p>
		<?php
		return;
	}
}

?>
<div class="commentHeader">
	<h4><?php _e('Comments', 'Arjuna'); ?></h4>
	<?php if(comments_open()): ?>
		<a href="#respond" class="btnReply btn"><span><?php _e('Leave a comment', 'Arjuna'); ?></span></a>
	<?php endif; ?>
	<?php if(pings_open()): ?>
		<a href="<?php trackback_url(); ?>" class="btnTrackback btn"><span><?php _e('Trackback', 'Arjuna'); ?></span></a>
	<?php endif; ?>
</div>
<?php if (have_comments()) { ?>
<ul class="commentList<?php if($arjunaOptions['commentDisplay']=='left') { echo ' commentListLeft'; } elseif($arjunaOptions['commentDisplay']=='right') { echo ' commentListRight'; } elseif($arjunaOptions['commentDisplay']=='alt') { echo ' commentListAlt'; } ?>">
	<?php wp_list_comments('callback=arjuna_get_comment'); ?>
</ul>

<?php arjuna_get_comment_pagination(); ?>

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
	<div class="replyHeader">
		<h4><?php _e('Leave a Comment', 'Arjuna'); ?></h4>
		<?php if (function_exists('cancel_comment_reply_link')) { ?>
			<div id="cancel-comment-reply" class="cancelReply"><?php arjuna_cancel_comment_reply_link('<span>'.__('Cancel Reply', 'Arjuna').'</span>');?></div>
		<?php } ?>
	</div>
	<?php if ( get_option('comment_registration') && !$user_ID ) : ?>
		<p style="margin-bottom:40px;"><?php printf(__('You must be %slogged in%s to post a comment.', 'Arjuna'), '<a href="'.get_option('siteurl').'/wp-login.php?redirect_to='.get_permalink().'">', '</a>'); ?></p>
	<?php else : ?>
		<form action="<?php echo get_option('siteurl'); ?>/wp-comments-post.php" name="reply" method="post" id="commentform">
			<input type="hidden" id="replyNameDefault" value="<?php _e('Your name', 'Arjuna'); ?>" />
			<input type="hidden" id="replyEmailDefault" value="<?php _e('Your email', 'Arjuna'); ?>" />
			<input type="hidden" id="replyURLDefault" value="<?php _e('Your website', 'Arjuna'); ?>" />
			<input type="hidden" id="replyMsgDefault" value="<?php _e('Your comment', 'Arjuna'); ?>" />
			<?php if ( $user_ID ): ?>
				<div class="replyLoggedIn">
				<?php printf(__('Logged in as %s.', 'Arjuna'), '<a href="'.get_option('siteurl').'/wp-admin/profile.php">'.$user_identity.'</a>'); ?> <a class="btnLogout btn" href="<?php echo wp_logout_url(get_permalink()); ?>" title="<?php _e('Log out of this account', 'Arjuna'); ?>"><span><?php _e('Logout', 'Arjuna'); ?></span></a>
				</div>
			<?php else : ?>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author)): ?> inputIA<?php endif; ?>" id="replyName" name="author" value="<?php if(!empty($comment_author)) { echo $comment_author; } else { _e('Your name', 'Arjuna'); } ?>" /></div>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author_email)): ?> inputIA<?php endif; ?>" id="replyEmail" name="email" value="<?php if(!empty($comment_author)) { echo $comment_author_email; } else { _e('Your email', 'Arjuna'); } ?>" /></div>
				<div class="replyRow"><input type="text" class="inputText<?php if(empty($comment_author_url)): ?> inputIA<?php endif; ?>" id="replyURL" name="url" value="<?php if(!empty($comment_author_url)) { echo $comment_author_url; } else { _e('Your website', 'Arjuna'); } ?>" /></div>
			<?php endif; ?>
			<?php
			if (function_exists('cancel_comment_reply_link')) { 
				//2.7 comment loop code
				comment_id_fields();
			}
			?>
			<input type="hidden" name="comment_post_ID" value="<?php echo $id; ?>" />
			<div class="replyRow"><textarea class="inputIA" id="replyMsg" name="comment"><?php _e('Your comment', 'Arjuna'); ?></textarea></div>
			<div class="replySubmitArea">
				<a href="<?php echo get_post_comments_feed_link(); ?>" class="btnSubscribe btn"><span><?php _e('Subscribe to comments', 'Arjuna'); ?></span></a>
				<button type="submit" class="inputBtn" value="Submit" name="submit"><?php _e('Leave comment', 'Arjuna'); ?></button>
			</div>
			<?php do_action('comment_form', $post->ID); ?>
		</form>
	<?php endif; // If registration required and not logged in ?>
</div>
<?php endif; ?>
