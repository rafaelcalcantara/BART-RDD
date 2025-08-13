res.list <- c("k1_1_k2_1_k3_0_k4_0.1_k5_0_p_2_rho_0.5",
              "k1_1_k2_1_k3_0_k4_0.1_k5_0_p_4_rho_0",
              "k1_1_k2_1_k3_0_k4_0.1_k5_1_p_2_rho_0",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_0_p_4_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_2_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_4_rho_0")
for (i in 1:length(res.list))
{
  dgp <- res.list[i]
  out <- matrix(0,100,5)
  for (sample in 1:100)
  {
    file <- paste0("sample_",sample,".csv")
    file.name <- paste0("Results/RMSE/",dgp,"/",file)
    out <- as.matrix(read.table(file.name))
    write.table(out[1,],paste0("Results/RMSE/",dgp,"/barddt_",file))
    write.table(out[2,],paste0("Results/RMSE/",dgp,"/tbart_",file))
    write.table(out[3,],paste0("Results/RMSE/",dgp,"/sbart_",file))
    write.table(out[4,],paste0("Results/RMSE/",dgp,"/polynomial_",file))
    write.table(out[5,],paste0("Results/RMSE/",dgp,"/ate_",file))
  }
}
