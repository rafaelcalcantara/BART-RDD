file.remove(dir(".")[grep("run",dir("."))]) ## Remove existing bash scripts for running sims
n <- 4000
s <- 10
k1 <- 2
k2 <- c(2,0.5)
k3 <- c(0,1)
k4 <- c(0.1,1)
k5 <- c(0.25,2)
p <- c(2,5)
rho <- c(0,0.7)
dgp <- matrix(0,14,7)
dgp <- matrix(1,7,7)+diag(1,7)
dgp <- rbind(dgp,matrix(2,7,7)+diag(-1,7))
dgp[,1] <- 1
for (i in 1:nrow(dgp))
{
  config <- dgp[i,]
  dgp[i,] <- c(k1[config[1]],k2[config[2]],k3[config[3]],k4[config[4]],
               k5[config[5]],p[config[7]],rho[config[6]])
}
models <- c("leaf.rdd","tbart","sbart","polynomial")
##
txt <- rep(NA,9)
txt[1] <- "#! /bin/bash"
txt[2] <- "#SBATCH --nodes=1                    # Run all processes on a single node"
txt[3] <- "#SBATCH --cpus-per-task=10            # Number of CPU cores per task"
txt[4] <- "#SBATCH --mail-type=END,FAIL # Mail events (NONE, BEGIN, END, FAIL, ALL)"
txt[5] <- "#SBATCH --mail-user=rafael.campellodealcantara@mccombs.utexas.edu # Where to send mail"
txt[6] <- "#SBATCH --mem=10G # Total memory limit"
txt[7] <- "#SBATCH --time=04:00:00 # Time limit hrs:min:sec"
txt[8] <- "#SBATCH --nice"
for (i in 1:nrow(dgp)) {
  txt[9] <- paste0("Rscript --verbose simulation_master.R ", paste(dgp[i,],collapse=" ")," ",
                   n," ",s," ",paste(models, collapse = " ")," > outputFile.Rout 2> errorFile.Rout", " &")
  writeLines(txt,paste0("run_sims",i,".sh"))
}