cutpoints <- (p+1)*sapply(1:s, function(i) sum(c-3*h[i] <= x[,i] & x[,i] <= c+3*h[i])) + 1
sapply(constraint.fails.1,summary)
sapply(constraint.fails.2,summary)
hist(constraint.fails.1[[sample]]/((ntrees_con+ntrees_mod)),
     main="Avg. no. of times condition i (Omin) was violated per tree\nfor each sample",
     xlab=paste0("Sample ",sample))
hist(constraint.fails.2[[sample]]/(ntrees_con+ntrees_mod),
     main="Avg. no. of times condition ii ('force split') was met per tree\nfor each sample",
     xlab=paste0("Sample ",sample))
cutoff.nodes.con
invalid.nodes.1.con
invalid.nodes.2.con
cutoff.nodes.mod
invalid.nodes.1.mod
colSums(sapply(invalid.nodes.2.con,rowSums))
colSums(sapply(invalid.nodes.2.mod,rowSums))
# hist(cutoff.nodes[[sample]]/(ntrees_con+ntrees_mod),
#      main="Nodes used for prediction per tree",
#      xlab=paste0("Sample ",sample))
# hist(invalid.nodes.2[[sample]]/cutoff.nodes[[sample]],
#      main="Prop. nodes used for prediction that break 'force split' condition",
#      xlab=paste0("Sample ",sample))