#!/bin/bash
nice Rscript --verbose generate_bash_scripts.R
scripts=$(find . -name "*run*")
for i in $scripts
do
chmod u+x $i
$i &
sleep 3
done
