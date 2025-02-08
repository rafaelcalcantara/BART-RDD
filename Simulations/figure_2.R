easy <- readRDS("easy_fits.rds")
hard <- readRDS("hard_fits.rds")
####
colors <- c("magenta3","dodgerblue", "green3","gold1")
colors[1] <- grDevices::adjustcolor(colors[1],alpha=1)

bgcolors <- colors

lincol <- adjustcolor('red',alpha=0.75)

colvec <- c(rep("black",3),"black",lincol)
cex = 1.75
pch <- c(rep(21,3),rep(21,2))
lwd <- c(rep(1,4),2)
pdf("fits.pdf",width=12,height=10)
layout(matrix(c(1,4,7,2,5,7,3,6,7),ncol=3),heights = c(2,2,1))
## Easy setup
mat <- easy
ylim <- c(min(mat[,-1])-0.1,max(mat[,-1])+0.1)
mat <- mat[order(mat[,1]),]
mat <- mat[,c(1,4,5,6,3,2)]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("p","n","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","p","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","n","p","p","l"),lty=1,lwd=lwd)
## Hard setup
mat <- hard
mat <- mat[order(mat[,1]),]
mat <- mat[,c(1,4,5,6,3,2)]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("p","n","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","p","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","n","p","p","l"),lty=1,lwd=lwd)
plot.new()
legend("center",legend=c("BARDDT","T-BART","S-BART","Polynomial"),
       col="black",pt.bg=colors[c(4,1,2,3)],pch=21,cex=1.25,title="Estimator",title.font = 2,ncol=4)
dev.off()