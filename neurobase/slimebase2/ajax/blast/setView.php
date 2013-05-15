<?php
session_start();
echo($_SESSION['view']);
$_SESSION['view'] = $_GET['view'];
echo($_SESSION['view']);
?>