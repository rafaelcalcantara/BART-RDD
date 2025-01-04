split <- 0.5
split1 <- 0.3
split2 <- 0.7
## Plots
par(mfrow=c(1,2))
plot(0:1,0:1,type="n",xlab="X",ylab="W",yaxt="n",xaxt="n",main="BARDDT")
xlim <- par("usr")[1:2]
ylim <- par("usr")[3:4]
axis(1,c(xlim[1],0.5,xlim[2]),labels = c("c-h","c","c+h"),tick = F)
segments(xlim[1],split,xlim[2],split)
px <- c(rep(xlim[1]+(0.5-xlim[1])/2,4),
        rep(0.5+(xlim[2]-0.5)/2,4))
py <- c(ylim[1]+(split1-ylim[1])/2,
        split1+(split-split1)/2,
        split+(split2-split)/2,
        split2+(ylim[2]-split2)/2)
py <- rep(py,2)
points(px,py,pch=as.character(1:8))
###
plot(0:1,0:1,type="n",xlab="X",ylab="W",yaxt="n",xaxt="n",main="T-BART")
axis(1,c(xlim[1],0.5,xlim[2]),labels = c("c-h","c","c+h"),tick = F)
abline(v=0.5)
segments(xlim[1],split1,0.5,split1)
segments(0.5,split2,xlim[2],split2)
points(px,py,pch=as.character(1:8))