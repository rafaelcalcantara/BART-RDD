setwd("~/Git/XBCF-RDD/")
## Setup
s      <- 1000
sample <- c(500,1000)
model  <- 3:6
xi <- nu <- kappa <- c(0.25,2)
## CGS
res.mat <- as.data.frame(matrix(0,32,5,dimnames=list(1:32,c("CGS","Model","Xi","Nu","Kappa"))))
rmse.cgs <- res.mat
cov.cgs <- res.mat
length.cgs <- res.mat
alpha.cgs <- res.mat
index <- 0
for (j in 1:length(model))
{
    for (k in 1:length(xi))
    {
        for (l in 1:length(nu))
        {
            for (m in 1:length(kappa))
            {
                index <- index+1
                dgp <- c(model[j],xi[k],nu[l],kappa[m])
                file <- paste0("Results/cgs_",500,"_",model[j],"_",xi[k],"_",nu[l],"_",kappa[m])
                cgs <- readRDS(paste0(file,".rds"))
                cgs <- lapply(cgs, function(x) x$pred)
                a <- sapply(cgs, function(x) mean(x<=xi[k]))
                alpha.cgs[index,] <- c(mean(a),dgp)
                cgs <- t(sapply(cgs,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
                rmse.cgs[index,] <- c(sqrt(mean((cgs[,1]-xi[k])^2)),dgp)
                cov.cgs[index,] <- c(mean(cgs[,2] <= xi[k] & xi[k] <= cgs[,3]),dgp)
                length.cgs[index,] <- c(mean(cgs[,3] - cgs[,2]),dgp)
            }
        }
    }
}
rm(cgs,res.mat)
save.image("Tables/cgs.RData")