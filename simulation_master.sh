#! /usr/bin/bash
for sig in 0.25 1 4
do
    for k in 0.1 0.5 1
    do
	echo "sig=$sig; k=$k"
	Rscript R/simulation_data.R $sig $k
	echo "Saved data for: sig=$sig; k=$k"
	echo "Running regressions: sig=$sig; k=$k"
	Rscript R/simulation_master.R
	echo "Saving simulation results: sig=$sig; k=$k"
	Rscript R/simulation_save_results.R $sig $k
	echo "Simulation finished: sig=$sig; k=$k"
    done
done
