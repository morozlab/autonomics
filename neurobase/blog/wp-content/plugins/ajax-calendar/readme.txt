=== Plugin Name ===
Contributors: johnny5
Donate link: http://urbangiraffe.com/about/support/
Tags: ajax, calendar, widget
Requires at least: 2.0
Tested up to: 2.9
Stable tag: trunk

AJAX Calendar is a plugin that will display an AJAXified WordPress calendar.

== Description ==

AJAX Calendar is a plugin that will display an AJAXified WordPress calendar. This enhances the functionality of the standard calendar by:

* Allowing the asynchronous navigation of months, without updating the page
* Adding a drop-down list of post titles in each month (through the «-» button)
* Adding full caching support to reduce server load

== Installation ==

The plugin is simple to install:

1. Download the zip file
1. Unpack the zip. You should have a directory called `ajax_calendar`, containing several PHP files
1. Upload the `ajax_calendar` directory to the `wp-content/plugins` directory on your WordPress installation. It is important you retain the `ajax_calendar` directory structure
1. Make any modifications to your theme, as required
1. Activate plugin

You can find full details of installing a plugin on the [plugin installation page](http://urbangiraffe.com/articles/how-to-install-a-wordpress-plugin/).

== Frequently Asked Questions ==

= Does this plugin provide any scheduling functionality? =

No, the calendar just shows your posting history.

== Screenshots ==

1. Example screenshot

== Documentation ==

Full documentation can be found on the [AJAX Calendar](http://urbangiraffe.com/plugins/ajax-calendar/) page.

== Changelog ==

= 2.1   =
* Swap from xajax to prototype.
* Fix cache deletion problems.

= 2.2   =
* Use wpurl instead of home

= 2.3   =
* Require WordPress 2.1+.
* Fix some problems.
* Use prototype lite

= 2.3.1 =
* Add option for always showing full details

= 2.4.0 =
* Use moo.tools for size
* use WP Object cache

= 2.4.1 =
* Added base CSS class
* Improved widget

= 2.4.2 =
* Use microajax for even better size and no clashes with prototype

= 2.4.3 =
* Fix validation errors

= 2.4.4 =
* Apply filters to titles

= 2.4.5 =
* Update widget class for 2.1 support

= 2.4.6 =
* WP 2.5 support

= 2.4.7 =
* WP 2.6

= 2.4.8 =
* Add support for category selection

= 2.4.9 =
* Fix DB prefix bug

= 2.4.10 =
* Fix category display