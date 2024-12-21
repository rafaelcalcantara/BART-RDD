setwd("~/Git/BART-RDD")
dgps <- 1:length(list.files("Data"))
txt <- rep(NA,9)
txt[1] <- "#! /bin/bash"
txt[2] <- "#SBATCH --nodes=1                    # Run all processes on a single node"
txt[3] <- "#SBATCH --cpus-per-task=10            # Number of CPU cores per task"
txt[4] <- "#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)"
txt[5] <- "#SBATCH --mail-user=rafael.campellodealcantara@mccombs.utexas.edu # Where to send mail"
txt[6] <- "#SBATCH --mem=10G # Total memory limit"
txt[7] <- "#SBATCH --time=04:00:00 # Time limit hrs:min:sec"
txt[8] <- "#SBATCH --nice"
for (i in 4:12) {
  txt[9] <- paste0("Rscript simulation_master.R ", i, " &")
  writeLines(txt,paste0("R/run_sims",i,".sh"))
}
