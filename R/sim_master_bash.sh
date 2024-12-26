#!/bin/bash
rm ../Data -r
rm ../Results -r
nice Rscript --verbose simulation_data.R
nice Rscript --verbose gen_bash_scripts.R
scripts=$(find . -name "*run*")
for i in $scripts
do
chmod u+x $i
$i &
sleep 3
done
