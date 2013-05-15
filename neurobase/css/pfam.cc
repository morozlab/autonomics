
/* pfam.css
   jt6 20060406 WTSI

   Stylesheet with generally applicable rules for the whole site.

   $Id: pfam.css,v 1.85 2009-10-29 14:40:18 jt6 Exp $
   
   Copyright (c) 2007: Genome Research Ltd.
   
   Authors: Rob Finn (rdf@sanger.ac.uk), John Tate (jt6@sanger.ac.uk)
   
   This is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with
   this program. If not, see <http://www.gnu.org/licenses/>.
    
*/   

/* global or widely used classes */

body {
  font-family: Verdana, Arial, Helvetica, sans-serif;
  font-size: small;
  color: #404040;
  background: #FFF;
  margin: 0;
  padding: 0;
}

.cleaner {
  clear: both;
  height: 1px;
  margin: 0;
  padding: 0;
}

ul, dl, ol {
  position: relative;
}

.ext {
  background: url(/static/images/external.gif) no-repeat center right;
  padding-right: 12px;
}

.feed {
  background: url(/shared/images/feed-icon-14x14.png) no-repeat center right;
  padding-right: 16px;
}

a:link,
a:visited,
a:focus,
a:active {
  color: #074987;
}

.hover,
a:hover {
  color: #666 ! important;
  text-decoration: underline ! important;
}

img {
  border: none;
}

img.info {
  vertical-align: middle;
}

.odd {
  background: #f3f3f3;
}

.even {
  background: #e7e7e7;
}

.small {
  font-size: smaller;
}

.nonbold {
  font-weight: normal;
}

.wide {
  width: 100%;
}

h1 {
  font-size: large;
  margin: 0.8em 0;
}
h2 {
  font-size: medium;
  margin: 0.5em 0;
}
h3 {
  font-size: small;
  margin: 0.4em 0;
}
h4, h5, h6 {
  font-size: x-small;
  margin: 0.2em 0;
}

/* hide the iframe that's used by the history manager in IE*/
#yui-history-iframe {
  position: absolute;
  top: 0;
  left: 0;
  width: 1px;
  height: 1px;
  visibility: hidden;
}

/* ========================================
 * header classes
 *
 * the header is laid out using "The Holy Grail" method from ALA, with
 * IE7 fixes from http://www.gerd-riesselmann.net/development/the-holy-grail-css-layout-fix-for-ie7
 */

#header {
  padding: 0 250px 0 250px;
  min-width: 400px;
}

#header .column {
  position: relative;
  float: left;
}

#navbar { /* centre */
  width: 100%;
}

div.siteLogoHolder { /* left */
  width: 250px;
  margin-left: -100%;
  left: 250px;
  padding: 10px 0 0 10px;
}

div.pfamLogoHolder { /* right */
  width: 250px;
  margin-right: -250px;
  padding: 10px 0 0 0;
}

#header > div.siteLogoHolder {
  left: -250px;
}

/* we need the right "column" to be the same width as the left in order to 
 * centre the navbar, but then we need to push the contents of the right column
 * back right to account for that */
div.pfamLogoHolder a,
div.pfamLogoHolder #searchBox {
  margin-left: 100px;
}

#siteLogo {
  background-image: url(/shared/images/sanger_logo_light.png);
/*  background-image: url(jfrc_gray_logo.gif); */
  background-repeat: no-repeat;
  width: 218px;
  height: 70px;
}

#pfamLogo {
  background: #FFF url(/static/images/pfam_logo.gif) no-repeat;
  width: 140px;
  height: 50px;
}

img.prePfamLogo {
  background-image: url(/static/images/pre_pfam_logo.gif) ! important;
}

img.h3PfamLogo {
  background-image: url(/static/images/h3_pfam_logo.gif) ! important;
}

/* tweak the site-wide search box so that it fits nicely under the Pfam logo */
#siteSearchField {
  width: 92px;
  margin-left: 3px;
}

/* navbar styles */
#navbar {
  font-size: 0.9em;
}

#navbar ul {
  margin: 40px 0 0 0;
  padding: 0;
  text-align: center;
}

#navbar li { 
  display: inline;
  letter-spacing: 2px;
}

#navbar span,
#navbar a {
  font-weight: bold;
}

#navbar a { 
  text-decoration: none;
  color: #105a94;
}

#navbar a.here { 
  color: #666;
}

/* a little warning message */

#jsWarning {
  border: 2px solid #600;
  background: #900;
  color: #FFF;
  margin: 0.4em auto;
  width: 50%;
  padding: 0.4em;
}

/* the "add annotation" button image */
img.addAnnotation {
	/* width: 100px; */
  height: 24px;
	vertical-align: middle;
  padding-bottom: 3px;
	border-style: none;
}

/* widgets in the jump box */

.entryField {
  border: 1px inset #074987;
  background: #cbdced url(/shared/images/search_gradient.png) repeat-x top left;
  margin: 1px 0;
  padding: 2px;
  color: #074987;
  font-weight: bold;
  font-size: x-small;
  float: left;
}
.disabledEntryField {
  background-color: #e7e7e7;
  color: #666;
}

/* jump box styles */

/* set the size of the field explicitly for the jumpBox */
#jumpBox .entryField {
  width: 110px; /* set to 145px for the right sidebar on the index page */
}

div.jumpError {
  clear: both;
  font-size: smaller;
  font-weight: bold;
  color: #600;
}

div.jumpSpinner {
  background-image: url(/shared/images/loading.gif);
  background-repeat: no-repeat;
  background-position: center left;
  padding-left: 16px;
  clear: both;
  font-size: smaller;
}

.glassButton {
  border: none;
  vertical-align: middle;
  margin: 1px;
  cursor: pointer;
  float: left;
  font-size: 1px;
  background-color: transparent;
  background-repeat: no-repeat;
  background-position: right;
}

.goButton {
  width: 33px;
  height: 20px;
  background-image: url(/static/images/go.png);
}
.disabledGoButton {
  width: 33px;
  height: 20px;
  background-image: url(/static/images/disabledGo.png);
}

.exampleButton {
  width: 66px;
  height: 20px;
  background-image: url(/static/images/exampleButton.png);  
}

/* ajax loading image */
.loading {
  background-image: url(/shared/images/loading.gif);
  background-repeat: no-repeat;
  background-position: center left;
  padding-left: 16px;
}

/* ---------------------------------------- */
/* search errors page */

.tips {
	margin: 0.2em auto;
  width: 60em;
}

.tips h2 {
	margin-top: 1em;
}

.tips dt {
  float: left;
  font-weight: bold;
	text-align: right;
  width: 15em;
}

.tips dd {
  margin-left: 16em;
	margin-bottom: 0.4em;
}

.tips p.back {
  text-align: center;
	font-size: large;
}

.tips p.back a {
	font-size: large;
}

/* ======================================== */
/* the footer. We use a wrapper div (#contentContainer) for the content
   when the page contains tabs, so add a separate style for that which 
   adds the background image, etc. */

#siteFooter {
  margin: 0.4em 0 0 0;
  padding: 0 0 0.4em 0;
  text-align: center;
  font-size: x-small;
}

#siteFooter a {
  text-decoration: none;
  font-weight: bold;
}

#siteFooter p {
  margin: 0.2em;
}

/* the footer for full page templates */
#footer span {
  margin: 0;
  padding: 0 0.4em 0 0;
}

#siteFooter p.spaced,
#footer p.spaced {
  letter-spacing: 3px;  
}

/* the footer for tab-layout templates */
#contentContainer #footer {
  clear: both;
  position: relative;
  width: 85%;
  float: right;
  border: none;
  border-top: 1px solid #074897;
  height: 1px;
}

/* ======================================== */
/* error page styles */

div.errorReport {
  width: 50%;
  margin: 1em auto;
  background: #EEE;
  border: 1px solid #CCC;
}

div.errorReport p {
  padding: 0.5em;
  margin: 0;
}

div.errorReport p.message {
  font-weight: bold;
  color: #833;
  margin-left: 2em;
}

/* ---------------------------------------- */
/* "details" table. Used in various places */

table.details {  
  margin: 0;
  padding: 0;
  border-spacing: 1px;
}

table.details thead tr {
  background: #BBB;
}

table.details th {
  padding: 0.2em;
/*  border: 2px solid #FFF; */
}

table.details td {
  vertical-align: middle;
  padding: 0.2em;
  text-align: center;
}

table.details td.label {
  width: 14em;
  padding: 0.2em;
  background: #CCC;
  text-align: right;
  font-weight: bold;
}

table.details tr.newRow td {
  margin-top: 2px;
  border-top: 1px solid #BBB;
}

td.stripeHover {
  background: #c2d2e1;
}

/* these are redundant...
 * table.details tr.odd {
 *   background: #CCC;
 * }
 * 
 * table.details tr.even {
 *   background: #DDD;
 * }
 */

table.details tr.summaryRow {
  background: #AAA;
}

form > div > input {
  margin: 0 1px;
}

/* these two rules are for the details table when it is split into
   two separate columns */
div.floatLeft {
  float: left;
  margin-right: 0.4em;
}

tr.continuation {
  font-size: smaller;
  background: none;
}

td.left {
  text-align: left !important;
}

/* fix the width of some columns */
table.details td.fixed {
  width: 10em;
}

/* different flavours of details table... stop me when it's getting ridiculous */
table.links td {
  width: 30em;
  text-align: left;
}

/* fix the width of a "notes" column */
td.notes {
  width: 40%;
}

/* ======================================== */
/* tool pages */

body.tool {
  padding: 80px 2px 2px 2px;
  background: #FFF url(/static/images/about_back.png) left top no-repeat; 
  color: #15538e;
}

.closer {
  background: url(/shared/images/close_grey.gif) left no-repeat;
  padding-left: 11px;
  margin: 0.2em;
}

body.tool h1 {
  margin: 0;
  padding: 0.4em 0;
}

/* ======================================== */
/* search results */

.sortIndicator th {
  background-image: url(/shared/images/up-down.gif);
}

.sortcol {
  cursor: pointer;
  padding-left: 16px ! important;
  padding-right: 16px ! important;
  background-repeat: no-repeat;
  background-position: right center;
}
.sortasc {
  background-image: url(/shared/images/up.gif);
}
.sortdesc {
  background-image: url(/shared/images/down.gif);
}
.nosort {
  cursor: default;
  background-image: none ! important;
}

/* this is the auto-generated link for sorting the column */
a.sortheader {
  text-decoration: none;
  display: block;
}

/* affects the key */

.key {
  padding: 4px;
  margin: 6px;
  border: 1px solid #CCC;
  background-color: #EEE;
  text-align: center;
}

.key h2 {
  color: #074987;
  text-align: left;
  margin: 0.2em auto;
  padding: 0;
  font-size: large;
  font-weight: bold;
}

.key p {
  margin: 0.4em 0;
  text-align: left;
  line-height: 1.4em;
}

/* other things in the key div: lists */

.key ul {
  text-align: left;
  margin-top: 0.2em;
  margin-bottom: 0.2em;
}

.key li {
  margin-top: 0.2em;
  margin-bottom: 0.2em;
}

.key dd {
  margin: 0 0 0 6em;
  text-align: left;
}

/* colour the progress bars according to job class */
div.pfamASearch .progressBar {
  background-color: #CBD;
}

div.pfamBSearch .progressBar {
  background-color: #CDB;
}

/* the top and bottom of the colour gradient in the functional similarity
   search results */
.key .gradTop {
  font-weight: bold;
  color: #FFF;
  background: #008000;
  padding: 0 0.2em;
}

.key .gradBot {
  font-weight: bold;
  color: #FFF;
  background: #C00000;  
  padding: 0 0.2em;
}

/* affects the key table */

#keyTable {
  margin-top: 4px;
  margin-bottom: 4px;
  width: 50%;
  margin-left: 25%;
  margin-right: 25%;
  border-spacing: 1px;
}

#keyTable thead {
  background: #AAA ! important;
  color: #FFF;
  font-weight: bold;
}

#keyTable tr.odd {
  background-color: #DDD;
}

#keyTable tr.even {
  background-color: #CCC;
}

/* affects the results div */

#results {
  margin: 5px;
}

/* affects the results table */

.resultTable {
  clear: both;
/*   border: 1px solid #000070; */
  margin: 0 auto;
  border-spacing: 1px;
}

.resultTable thead tr {
   background: #074987;
   color: #FFF;
   font-weight: bold;
}
.resultTable thead a {
   color: #FFF;
   font-weight: bold;
}

.resultTable .rowNum {
  display: none;
}

.resultTable p {
  display: none;
}

.resultTable .fixedWidth {
  width: 6em;
}

.resultTable .desc {
  text-align: left;
}

.resultTable th {
  padding: 3px;
}

.resultTable td {
  padding: 3px;
  text-align: center;
  background-position: center center;
  background-repeat: no-repeat;
  background-image: url(/shared/images/blank_tick.gif);
}

.resultTable td.ticked {
  background-image: url(/shared/images/tick.gif);
}

.highlight {
  background: #c2d2e1 ! important; /* #CBD */
}

/* ======================================== */
/* rules that are common to several blocks */

div.domain {
  padding: 6px;
}

div.pdbImageFragment {
  float: right;
  margin: 0.2em;
  padding: 0.3em;
  width: 200px;
  font-size: 0.8em;
  overflow: hidden;
  border: 1px solid #074987;
  background: #e7edf3;
}

div.pdbImageFragment h1 {
  font-size: 1em;
  text-decoration: underline;
}

img.pdbImage {
  width: 200px;
  height: 200px;
}

.pdbTitle  {
  display: inline;
  color: #074987;
  font-weight: bold;
}

hr {
  background-color: #e7e7e7;
  border: none;
  height: 1px;
  width: 60%;
  margin: 2em auto;
}

hr.short {
  display: none;
}

.plainSequence {
  font-family: monospace;
  font-size: 0.8em;
  margin: 0.2em;
  padding: 0.2em;
  border: 1px solid #CCC;
  background: #DDD;
  float: left;
  clear: both;
}

.centreWrapper {
  float: left;
  position: relative;
  left: 50%;
}

.centredBlock {
  position: relative;
  left: -50%;
}

/* ---------------------------------------- */
/* pfam graphics */

div.pgholder {
  overflow: auto;
  /* position: relative; */
}

div.graphicRow {
  clear: both;
  padding: 0.2em;
}

span.graphicLabel {
  width: 24em;
  margin: 0;
  padding: 0;
}

div.graphicRow p,
div.graphicRow h3 {
  margin: 0.2em 0;
  padding: 0;
}

img.graphicImage {
  margin: 0.2em 0;
  padding: 0;
  float: left;
}

div.allSequences {
  overflow: auto;
  border: 1px solid #AAA;
  margin: 0.6em;
  max-height: 250px;
  position: relative;
}

div.allSequences .odd {
  background: #DED;
}

div.allSequences .even {
  background: #CDC;
}

#selectedAccs {
  width: 90%;
  margin: 0 auto;
}

#selectedAccs a,
#selectedAccs span {
  width: 6em;
  float: left;
  margin-bottom: 0.4em;
}

/* ---------------------------------------- */
/* a list of items */

span.listItem {
  float: left;
  padding-left: 0.8em;
  width: 12em;
}

/* a padded heading, which we use after the list of items... */
h2.padded {
  margin: 0;
  padding: 1em 0;
  clear: left;
}

/* ---------------------------------------- */
/* mimic a link using a span... */

.link {
  cursor: pointer;
  text-decoration: underline;
  color: #074987;
}

/* a slightly heavier link class... */
.moreLink {
  cursor: pointer;
  text-decoration: none;
  color: #074987;
  font-weight: bold ! important;
  text-decoration: underline ! important;
}

/* ---------------------------------------- */
/* tooltip styles */

.prototip .tooltip {
  background: #efefef;
  border: 1px solid #074987;
}

.prototip .tooltip .title {
  background: #074987;
  color: #cbdced;
  font-weight: bold;
  padding: 0 0.2em;
}

.prototip .tooltip .content {
  color: #074987;
  margin: 0.2em;
}

/* for the benefit of IE, otherwise it leaves a great wide gap
 * at the left of the <dl> */
.prototip dl {
  margin: 0;
}

.prototip p {
  margin: 0;
  padding: 0.2em;
}

.prototip div.scale {
  font-family: monospace;
  font-size: smaller;
  text-align: left;
}

.prototip div.scale span {
  margin: 0;
  padding: 0;
}

.prototip dt {
  padding: 0;
  font-weight: bold;
  float: left;
}

.prototip dd {
  margin-left: 8em;
  padding: 0 0.2em;
}

.prototip .inactive {
  color: #aaaaaa;
}

/* we know that the terms are very small in the tooltip on the alignments */
#alignmentKey dt {
  width: 6em ! important;
}

#alignmentKey dd {
  margin-left: 6em ! important;
}

/* ---------------------------------------- */
/* format the mirrors list and citation */

/* the citation on the index page */
#citation {
  width: 60%;
  float: left;
}

/* all other citations */

div.citation {
  padding-bottom: 0.8em;
}

div.citation span.title {
  font-style: italic;
}
div.citation span.ref {
  display: block;  
  padding: 0.2em 0;
}
div.citation span.jrnl {
  font-weight: bold;
}

/* mirrors list */

ul.mirrors {
  margin: 0;
  padding: 0;
  list-style: none;
}

ul.mirrors li {
  padding: 0.2em 0 ! important;
}

.flag {
  padding-left: 14px;
}

.mirror {
  padding-left: 14px;
  background-position: left;
  background-repeat: no-repeat;
}

/* the national flags */
.uk {
  background-image: url(/static/images/flags/british_flag.gif);
}
.us {
  background-image: url(/static/images/flags/american_flag.gif);
}
.se {
  background-image: url(/static/images/flags/swedish_flag.gif);
}
.kr {
  background-image: url(/static/images/flags/south_korea.gif);
}
.fr {
  background-image: url(/static/images/flags/french_flag.gif);
}

/* ======================================== */
/* a few rules for a page that apologises when the database is down */

#dbd {
  margin: 120px auto 0 auto;
  width: 485px;
  color: #074987;
  text-align: justify;
  font-size: 12px;
}

#dbd h1 {
  font-size: 20px;
  text-align: center;
}

#dbdSiteLogo {
  background: url(/shared/images/sanger_logo_light.png) no-repeat;
  width: 200px;
  height: 70px;
  left: 0;
  position: relative;  
}

#dbdPfamLogo {
  background: url(/static/images/pfam_logo_large.gif) no-repeat;
  width: 200px;
  height: 70px;
  right: 0;
  position: relative;  
}

/* ======================================== */
/* species tree styles */

.highlightSeed {
  background: #c2d2e1; /* #CBD */
}

td.nodeSummaryCell {
  width: 16em;
}

div.nodeSummary {
  margin-left: 0.2em;
}

div.nodeSummary div {
  float: left;
  border-width: 1px;
  border-style: solid;
  width: 3em;
  height: 1.2em;
  padding: 0.1em;
  margin: 0.1em;
  text-align: center;
}

.specSum {
  background: #FDD;
  border-color: #B99;
}
.seqSum {
  background: #DFD;
  border-color: #9B9;
}
.domSum {
  background: #DDF;
  border-color: #99B;
}

#treeDiv {
  min-height: 31em;
}

#treeTools {
  right: 3.2em;
  width: 18em;
  padding: 0.4em;
  border: 1px solid #AAA;
  background: #e8e8e8;
}

#speciesTreeDesc {
  margin-right: 20em;
}

div.blockContent #treeTools {
  position: fixed;
  cursor: move;
}

#treeTools ul {
  padding-left: 2em;
}

#treeTools ul.bare {
  list-style: none;
  margin: 0;
  padding: 0;
}

#treeTools li {
  padding: 0.2em 0;
}

#treeTools h3 {
  font-size: 1em;
  padding: 0.4em 0 0.2em 0;
  margin: 0;
}

#toolsHeader h3 {
  margin: 0;
  padding: 0.2em 0;
  float: left;
}

#toolsToggle {
  cursor: pointer;
  float: right;
  padding: 0.2em 0;
}

#stError {
  color: #F00000;
}

/* ======================================== */
/* about page */

#about {
  width: 700px;
  margin: 0 auto;
}

#about p {
  line-height: 1.6em;
}


/* ======================================== */