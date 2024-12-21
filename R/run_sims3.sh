#! /bin/bash
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --cpus-per-task=10            # Number of CPU cores per task
#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rafael.campellodealcantara@mccombs.utexas.edu # Where to send mail
#SBATCH --mem=10G # Total memory limit
#SBATCH --time=04:00:00 # Time limit hrs:min:sec
#SBATCH --nice
Rscript simulation_master.R 3 &
