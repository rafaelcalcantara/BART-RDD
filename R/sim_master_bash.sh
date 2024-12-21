#!/bin/bash
nice Rscript --verbose simulation_data.R
nice Rscript --verbose gen_bash_scripts.R
scripts=$(find . -name "*run*")
for i in $scripts
do
chmod u+x $i
$i &
sleep 1
done
