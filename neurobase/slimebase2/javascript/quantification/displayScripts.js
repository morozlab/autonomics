// JavaScript Document
function showPasteBox(){
	var pasteBox = document.getElementById('pasteBox');
	var fileBox = document.getElementById('fileBox');
	pasteBox.style.display = 'block';
	fileBox.style.display = 'none';
}

function showFileBox(){
	var pasteBox = document.getElementById('pasteBox');
	var fileBox = document.getElementById('fileBox');
	pasteBox.style.display = 'none';
	fileBox.style.display = 'block';	
}