#!/bin/csh
if ($#argv != 2) then
  echo "takes 2 args: file and command   expects file to be a list of files"
  exit
endif

foreach line ( "`cat $1`" )
#  echo "\n"
#  echo $2 $line
  $2 $line
end
