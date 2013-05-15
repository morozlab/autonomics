This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


---------------------------------------------------------------------


Please respect that we have put a healthy amount of time in the development, testing and the maintenance of this theme. You can use this theme free of charge according to the GPL license. However, we ask a small, yet completely optional, favor of you: If you like the theme, we would appreciate if you keep the credits in the footer of the theme intact. We have intentionally made them small and not intrusive. It's a small favor for a project that, had you hired someone to design and develop this theme for you, cost several thousand dollars. 


Translation
---------------------------------------------------


Help and how to translate Arjuna into your own language can be found here: http://192.168.2.100/srssolutions.com/en/downloads/arjuna_wordpress_theme#faq

Currently available translations:

- Italian (thanks to Gianni Diurno, gidibao.net)
- Lithuanian (thanks to Vaidas Juknevicius aka SeniZ)
- Chinese (thanks to Liuyue, liuyue.asia)
- French (thanks to Laurent Measson and Sebastien Violet)
- Brazilian (thanks to Pedro Spoladore)
- Spanish (thanks to José Marín)
- German (thanks to Frank Weichbrodt and Markus Liebl)
- Czech (thanks to Jirka Knapek and Ivan)
- Turkish (thanks to Serhat Yolaçan)
- Hungarian (thanks to Márton Bakos)
- Slovak


More Information and Support
---------------------------------------------------


Changelog: http://www.srssolutions.com/en/downloads/arjuna_wordpress_theme#changelog
FAQ: http://www.srssolutions.com/en/downloads/arjuna_wordpress_theme#faq
Report Bug: http://www.srssolutions.com/en/downloads/bug_report


Important notes about this theme:
---------------------------------------------------

* The two menus that are supported by this theme do NOT extend to more than one line. If you have more menu items and they don't fit, the theme has to be customized to accommodate your needs. Generally speaking, it is good practice to keep the options that a user of your website has limited. If you find that you have too many top categories, we advise to minimize the amount of top categories and create your category logic more intelligently and *user-friendly*. Besides, if you use sub categories, they will be shown in a dropdown menu as supported by Arjuna.


About IE6 Support:
---------------------------------------------------

Obviously, the IE6 version of this theme looks a little different than the version displayed in modern browsers. We had tried to simulate the modern version as much as possible but due to the nature of the IE6 rendering engine, we had to make certain decisions to cut down on some visual elements (layout remains completely the same!). The only alternative that we considered for a short time was to use GIF transparency, however, the graphics looked jagged at many areas so we decided to create an IE6 only version of the theme with a simpler background than the original. We believe the result is more than acceptable. Remember the horribly designed websites that had been created back in the days where IE6 was launched? This is definitely not one of them.
Nonetheless, we still hope that IE6 will rest in peace soon.


Changelog:
---------------------------------------------------

1.3.4
- Fix: The default setup with the sidebar having a normal width and being on the right caused issues with the display.
1.3.3
- Localization: German has been updated to 1.3.
- Localization: Brazilian has been updated to 1.3.2.
- Change: Improved the display of numbered and bulleted lists. Nested lists now show different list styles to help separate them visually. In addition, the margins have been altered slightly.
1.3.2
- Localization: Italian has been updated to 1.3.
- Fix: The new exclude pages/categories feature had some errors.
1.3
- New: It's now possible to hide comments completely, including the Comments section and the buttons, if comments are disabled for posts or pages.
- New: You can now exclude pages and categories from both header menus. It's pretty basic right now due a number of limitations imposed by the default WordPress functions. Therefore, the admin currently doesn't show any hierarchy information, but it should suffice for most purposes for now.
- New: You can now disable the info bar of static pages. This will render the title of static WordPress pages without the author name and the date it has been published. (Note: The info bar is turned off by default now!)
- New: Added a Facebook button to the sidebar button collection.
- New: You can now disable the labels next to the sidebar buttons (for the RSS, Twitter and Facebook button).
- New: Added two new header colors, Sea Green and Khaki.
- Change: The comment pagination has finally been implemented properly. It does not display any longer if no comment pages are available or if comment pagination has been disabled.
- Change: Added a few additional checks for the twitter URL. 
- Change: Changed the background color of the pagination bar on top of the posts to fit slightly better than the previous grayish background.
- Change: The names of the widget bars of the sidebar have been changed. Note: This might mean that you need to reassign your widgets if you are upgrading Arjuna.
- Fix: Bullet lists and numbered lists in the content of posts rendered with no white space below them.
- Fix: Some international blogs with a language that has characters not part of the standard ASCII table seemed to cause major rendering issues in several browsers, mostly Internet Explorer. Arjuna is now encoded in UTF-8.
- Localization: Slovak (1.2) has been added.
1.2.5
- Localization: Hungarian (1.2) has been added.
- Fix: The fix for the lightbox in 1.2.4 has caused the previous bug in 1.1.3, where a dropdown in the first header menu would appear behind the second header menu, to reappear again.
1.2.4
- Localization: Turkish (1.2) has been added.
- Localization: Italian has been updated to 1.2.3.
- Fix: The header does not appear above any Lightbox 2 layers anymore.
- Fix: A weird space character at the very beginning of the header.php file caused a bunch of issues, especially with Internet Explorer.
1.2.3
- Fix: The Edit in Admin button still appeared on pages even if you are not logged in.
- Fix: When IE6 optimization was turned on, the new navigation links within single posts did not render properly.
- Fix: The Edit in Admin button in IE6 was messed up.
1.2.1
- Fix: The Edit in Admin button always appeared even if nobody was logged in. It now only appears if you are logged in into the WordPress admin.
- Fix: If you chose to display pages in the second header menu, the items could not be sorted in descending order. Instead Arjuna always switched back to ascending.
- Fix: Some minor bugs in the admin.
- Localizations: German and French have been updated.
1.2
- New: Native support for pagination has been added.
- New: It's now possible to enable links to previous and next posts on permalink pages, i.e. the URL where one single post/page is displayed.
- New: The RSS button in the sidebar can now be disabled.
- New: Arjuna now integrates a simple Twitter icon, which will appear right next to the RSS icon in the sidebar.
- Localization: Czech (1.2) has been added.
- Localization: Lithuanian, Spanish, French, Brazilian, and German have been updated to 1.2.
- Change: The custom CSS now also works even if Arjuna has no write permissions to the theme directory. The CSS rules will be included in the HTML HEAD. This ensures maximum compatibility with a variety of setups.
- Change: Some minor performance optimizations.
- Fix: When the display of time in posts is disabled in the admin, the date would inaccurately append the time without any space or words in between.
- Fix: When the date format of comments was set to the default date format (showing a date instead of the elapsed time), the date was not correctly translated by WordPress. This only occured in non-English WordPress installations.
- Fix: When users had to be logged in to post a comment, the sidebar would be placed below the whole post.
- Fix: There had been some IE6 bugs since version 1.1.3, that appeared with and without IE6 optimization enabled.
1.1.4
- Localization: Brazilian has been updated to 1.1.3.
- Fix: When displaying the elapsed time of a comment, there was an error with handling time zones. That way, some comments showed a negative number of seconds.
- Fix: The feature to add custom CSS is now only enabled if the application has sufficient write permissions to create a user-style.css file.
1.1.3
- New: It's now possible to add your custom CSS rules via the admin. This will ensure that your custom CSS is not overwritten when you update Arjuna automatically.
- New: Posts and pages can now show not only the date but also the time of when the post/page has been published.
- New: Included support for H1 tags in posts.
- Localization: Brazilian has been added.
- Localization: German has been added.
- Localization: French has been updated.
- Localization: Spanish has been updated.
- Change: Optimized Arjuna for SEO purposes. The titles of single pages/posts are now H1 tags, while the titles on archive pages remain in H2 elements.
- Change: Updated .POT localization file.
- Fix: When someone submitted a comment and did not enter a website URL, the comment author would link to http://Yourwebsite.
- Fix: Full-width page templates still did not work properly on some setups.
- Fix: If an item in the first header menu had a dropdown of more than four items, then the dropdown would get cut at the bottom.
- Fix: Modified some minor things in header.php to make it W3C compliant.
1.1.2
- Localization: French has been added.
- Localization: Italian has been updated.
- Change: Designed the Edit link to look similar to the other buttons.
- Change: Updated .POT file.
- Fix: When using Arjuna with WPML, the Home link would always link to the default language (the WordPress root). It now displays the localized home.
- Fix: When some IE8 browsers emulate IE6 for compatibility reasons and Arjuna was set to use IE6 optimization, the website's layout would have some issues. 
1.1.1
- New: It's now possible to disable default widgets in the sidebar when the widget bars are empty. This allows for more flexibility in choosing which sidebars to use. For example, you can now disable default widgets and exclusively use the two-column or single-column sidebar, if you wish so.
- New: Added default styles for headings in posts (heading 2 to 4).
- Change: Improved the backend and the way options are retrieved.
- Fix: When some IE7 browsers emulate IE6 for compatibility reasons and Arjuna was set to use IE6 optimization, the website's layout would have some issues.
- Fix: Full-width pages did not work any longer for layouts other than the Arjuna default.
1.1
- New: The first header menu can now be extended over two rows of links.
- New: The buttons in the second header menu can now be visually separated.
- New: You can now select between two different header images. More will become available in the future. NOTE: As of this version, the layout and structure of the header has undergone significant changes.
- New: The Home button in the second header bar can now be disabled. In addition, the Home button shows an icon now.
- Localization: Chinese has been added.
- Localization: Italian has been updated.
- Change: Optimized JavaScript and reduced its file size.
- Fix: There were some non-localized strings in functions.php as well as some that PoEdit could not parse properly.
- Fix: If pagination is provided by the plugin "WP-paginate", it did not work correctly on archive pages.
- Fix: Bold text, italic text and lists (ordered and unordered) within a post did not display correctly.
- Fix: A two-level dropdown menu in a left-aligned first header menu was displaying incorrectly.
- Fix: IE6 displayed a considerably larger font size for widget bars that used the default font size, e.g. the WP calendar.
1.0.5
- New: It's now possible to display the author of a post in the post header.
- Change: Changed the display of dates to use the default date format specified in Settings > General instead of using a theme-based, localized date format. This allows for grammatically correct date formats in languages other than English where there is no translation available for Arjuna yet.
- Localization: Italian has been added.
- Localization: Lithuanian has been added.
1.0.4
- Fix: The bug from version 1.0.2 introduced another bug. We reverted the bug. Please refer to the note in the README.txt file.
1.0.3
- Fix: There were some issues with the sidebars when widget bars are deactivated.
- Fix: When the first header menu had been disabled, the space between the header and the content would have been twice as large.
- Fix: In the admin, when you selected custom for "Append to page title" and left the field empty, it would use the old value again. Arjuna now does not append anything.
- Fix: In admin, disabling the IE6 optimization feature was not possible.
- Fix: When the first header menu was displayed in IE6, the second header menu was slightly shifted upwards.
- Fix: The bottom of the sidebar in IE6 rendered a small part of the RSS icon.
- Fix: If the content of a widget bar in the sidebar explicitly extended over the width of the sidebar, e.g. in the case of an image, the sidebar in IE6 would break completely.
1.0.2
- Fix: There was a small issue with the header bars where they would overrun the header area if too many items were shown.