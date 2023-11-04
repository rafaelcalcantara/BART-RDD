## Reading command line arguments to save results
args <- commandArgs(trailingOnly=TRUE)
sig <- as.numeric(args[1])
k <- as.numeric(args[2])
##
bart.rdd.1 <- readRDS("Results/bart_rdd_1.R")
bart.rdd.2 <- readRDS("Results/bart_rdd_2.R")
bart.rdd.3 <- readRDS("Results/bart_rdd_3.R")
bcf.1 <- readRDS("Results/bcf_1.R")
bcf.2 <- readRDS("Results/bcf_2.R")
bcf.3 <- readRDS("Results/bcf_3.R")
sbart.1 <- readRDS("Results/sbart_1.R")
sbart.2 <- readRDS("Results/sbart_2.R")
sbart.3 <- readRDS("Results/sbart_3.R")
tbart.1 <- readRDS("Results/tbart_1.R")
tbart.2 <- readRDS("Results/tbart_2.R")
tbart.3 <- readRDS("Results/tbart_3.R")
ckt.1 <- readRDS("Results/ckt_1.R")
ckt.2 <- readRDS("Results/ckt_2.R")
ckt.3 <- readRDS("Results/ckt_3.R")
res <- list(bart.rdd.1=bart.rdd.1,bart.rdd.2=bart.rdd.2,bart.rdd.3=bart.rdd.3,
            bcf.1=bcf.1,bcf.2=bcf.2,bcf.3=bcf.3,
            sbart.1=sbart.1,sbart.2=sbart.2,sbart.3=sbart.3,
            tbart.1=tbart.1,tbart.2=tbart.2,tbart.3=tbart.3,
            ckt.1=ckt.1,ckt.2=ckt.2,ckt.3=ckt.3)
saveRDS(res,paste0("Results/results_",sig,"_",k,".rds"))
