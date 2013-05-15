<?php

$nt = $_GET['sequence'];

$array = array();
exec("perl /var/www/nb_scripts/slimebase2translate.pl " . $nt, $array);

echo(json_encode($array));

?>