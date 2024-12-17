  set.seed(0)
  library(XBART)
  N <- c(500,1000,2500,5000)
  burnin <- 50
  num_sweeps <- 150
  max_depth <- 250
  Ow <- c(0.057,0.016,0.012,0.009)
  out <- list(`500`=NA,`1000`=NA,`2500`=NA,`5000`=NA)
  for (n in N)
  {
    x <- 2*rbeta(n,2,4)-0.75
    w <- runif(n,0,1)
    z <- x > 0
    y <- 3*x^5 - 2*x^2 + x + w^2 + z*(1 + w - mean(w) + (exp(x)-1)/(1+exp(x))) + rnorm(n,0,1)
    train <- -0.3 < x & x < 0.3
    x <- x[train]
    w <- w[train]
    z <- z[train]
    y <- y[train]
    if (n==500) Owidth <- Ow[1]
    if (n==1000) Owidth <- Ow[2]
    if (n==2500) Owidth <- Ow[3]
    if (n==5000) Owidth <- Ow[4]
    bart.parallel <- system.time({
      XBART(y,cbind(x,w,z),num_trees = 10, num_sweeps = num_sweeps,
            num_cutpoints = n, burnin = burnin, parallel = T, nthread = 11,
            max_depth = max_depth)
    })
    bart <- system.time({
      XBART(y,cbind(x,w,z),num_trees = 10, num_sweeps = num_sweeps,
            num_cutpoints = n, burnin = burnin, parallel = F,
            max_depth = max_depth)
    })
    bart.rdd.parallel <- system.time({
      XBCF.rd(y=y,W=w,X=x,c=0,
              Owidth = Owidth, Omin = 1, Opct = 0.95,
              num_cutpoints = n,
              num_sweeps = num_sweeps,
              num_trees_con = 5,
              num_trees_mod = 5,
              burnin = burnin,
              p_categorical_con = 0,
              p_categorical_mod = 0,
              parallel = T, nthread = 11,
              max_depth = max_depth)
    })
    bart.rdd <- system.time({
      XBCF.rd(y=y,W=w,X=x,c=0,
              Owidth = Owidth, Omin = 1, Opct = 0.95,
              num_cutpoints = n,
              num_sweeps = num_sweeps,
              num_trees_con = 5,
              num_trees_mod = 5,
              burnin = burnin,
              p_categorical_con = 0,
              p_categorical_mod = 0,
              max_depth = max_depth)
    })
    out[[as.character(n)]] <- rbind(bart.parallel,bart,bart.rdd.parallel,bart.rdd)
  }
  ##
  pmat <- t(sapply(out, function(mat) mat[,3]))
  matplot(pmat,type="b",pch=19,lty=1,ylab="Time (seconds)",xlab="N",xaxt="n",bty="n")
  axis(1,at=1:4,labels=rownames(pmat))
  legend("topleft",legend = colnames(pmat),col=1:4,lty=1,pch=19,ncol=2)
  matplot(pmat[,c(2,4)],type="b",pch=19,lty=1,ylab="Time (seconds)",xlab="N",xaxt="n",bty="n")
  axis(1,at=1:4,labels=rownames(pmat))
  legend("topleft",legend = colnames(pmat)[c(2,4)],col=1:2,lty=1,pch=19,ncol=2)