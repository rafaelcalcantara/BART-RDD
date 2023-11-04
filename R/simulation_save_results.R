## Reading command line arguments to save results
args <- commandArgs(trailingOnly=TRUE)
sig <- as.numeric(args[1])
k <- as.numeric(args[2])
##
bart.rdd.1 <- readRDS("Results/bart_rdd_1.rds")
bart.rdd.2 <- readRDS("Results/bart_rdd_2.rds")
bart.rdd.3 <- readRDS("Results/bart_rdd_3.rds")
bcf.1 <- readRDS("Results/bcf_1.rds")
bcf.2 <- readRDS("Results/bcf_2.rds")
bcf.3 <- readRDS("Results/bcf_3.rds")
sbart.1 <- readRDS("Results/sbart_1.rds")
sbart.2 <- readRDS("Results/sbart_2.rds")
sbart.3 <- readRDS("Results/sbart_3.rds")
tbart.1 <- readRDS("Results/tbart_1.rds")
tbart.2 <- readRDS("Results/tbart_2.rds")
tbart.3 <- readRDS("Results/tbart_3.rds")
ckt.1 <- readRDS("Results/ckt_1.rds")
ckt.2 <- readRDS("Results/ckt_2.rds")
ckt.3 <- readRDS("Results/ckt_3.rds")
res <- list(bart.rdd.1=bart.rdd.1,bart.rdd.2=bart.rdd.2,bart.rdd.3=bart.rdd.3,
            bcf.1=bcf.1,bcf.2=bcf.2,bcf.3=bcf.3,
            sbart.1=sbart.1,sbart.2=sbart.2,sbart.3=sbart.3,
            tbart.1=tbart.1,tbart.2=tbart.2,tbart.3=tbart.3,
            ckt.1=ckt.1,ckt.2=ckt.2,ckt.3=ckt.3)
saveRDS(res,paste0("Results/results_",sig,"_",k,".rds"))
