### ADD STRAIGHT BIAS
setwd("~/Git/BART-RDD/")
load("Results/sims.RData")
## Tables
stargazer::stargazer(cbind(params[,-4],round(rmse,3)),summary=F,rownames = F)
# stargazer::stargazer(cbind(params[,-4],round(cov,3)),summary=F,rownames = F)
# stargazer::stargazer(cbind(params[,-4],round(ci,3)),summary=F,rownames = F)
stargazer::stargazer(cbind(params[,-4],round(rmse.cate,3)),summary=F,rownames = F)
stargazer::stargazer(cbind(params[,-4],round(cov.cate,3)),summary=F,rownames = F)
stargazer::stargazer(cbind(params[,-4],round(ci.cate,3)),summary=F,rownames = F)
cov.tab <- cbind(params[,-4],round(cov,3))
ci.tab <- cbind(matrix(NA,nrow(ci),3),round(ci,3))
ci.tab[,-c(1,2,3)] <- apply(ci.tab[,-c(1,2,3)],2,function(i) paste("(",i,")",sep=""))
cov.ci.tab <- do.call("rbind",lapply(1:nrow(cov), function(i) rbind(cov.tab[i,],ci.tab[i,])))
stargazer::stargazer(apply(cov.ci.tab,2,as.character),summary=F,rownames=F,
                     label="tab:cov.ci",
                     title="Coverage rate and interval sizes (in parenthesis) for the ATE",
                     font.size = "small",column.sep.width = "0pt",
                     out="C:/Users/Rafael Alcantara/Dropbox/RDD-BART/Tables/sim_cov_ci_ate.tex")
## Plots
par(bty="L")
### Overall
ticks <- paste0("(",do.call(paste,c(params[,1:3],sep=",")),")")
matplot(rmse,lty=1,pch=19,type="b",ylab="RMSE",xaxt="n",
        cex.axis=0.75,cex.lab=0.75,ylim=c(min(rmse),max(rmse)+0.1),
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.65)
legend("topleft",legend=colnames(rmse),col=1:ncol(rmse),bty="n",lty=1,pch=19,
       ncol=3,cex=15)
####
matplot(abs(sb),lty=1,pch=19,type="b",ylab="Absolute bias",xaxt="n",
        cex.axis=0.75,cex.lab=0.75,ylim=c(min(abs(sb)),max(abs(sb))+0.1),
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.65)
legend("topleft",legend=colnames(sb),col=1:ncol(v),bty="n",lty=1,pch=19,
       ncol=3,cex=15)
####
matplot(v,lty=1,pch=19,type="b",ylab="Variance",xaxt="n",
        cex.axis=0.75,cex.lab=0.75,ylim=c(min(v),max(v)+0.1),
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.65)
legend("topleft",legend=colnames(v),col=1:ncol(v),bty="n",lty=1,pch=19,
       ncol=3,cex=15)
####
plot(colMeans(ci),colMeans(cov),col=1:6,pch=19,bty="L",
     xlab="Average interval length",ylab="Average coverage rate",
     cex.axis=0.75,cex.lab=0.75)
legend("topleft",legend=colnames(cov),title="Estimator",ncol=2,
       cex=0.55,pt.cex=0.75,title.cex=15,
       pch=19,col=1:6)
### RMSE
rmse.plot <- reshape(cbind(params,rmse),direction="long",varying=list(5:9),timevar="Model",times=colnames(rmse),v.names = "RMSE")
rmse.plot$Model <- factor(rmse.plot$Model, levels=colnames(rmse))
#### delta_mu
pdf("Figures/rmse_delta_mu.pdf")
par(bty="L")
boxplot(rmse.plot$RMSE~rmse.plot$delta.mu+rmse.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="RMSE",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(rmse)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(rmse)-1),par("usr")[3]-0.03, labels = colnames(rmse),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(rmse)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.mu),title=expression(delta[mu]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### delta_tau
pdf("Figures/rmse_delta_tau.pdf")
par(bty="L")
boxplot(rmse.plot$RMSE~rmse.plot$delta.tau+rmse.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="RMSE",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(rmse)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(rmse)-1),par("usr")[3]-0.03, labels = colnames(rmse),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(rmse)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.tau),title=expression(delta[tau]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### tau
pdf("Figures/rmse_tau.pdf")
par(bty="L")
boxplot(rmse.plot$RMSE~rmse.plot$tau+rmse.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="RMSE",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(rmse)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(rmse)-1),par("usr")[3]-0.03, labels = colnames(rmse),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(rmse)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$tau),title=expression(tau),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### Bias
sb.plot <- reshape(cbind(params,sb),direction="long",varying=list(5:9),timevar="Model",times=colnames(sb),v.names = "sb")
sb.plot$Model <- factor(sb.plot$Model, levels=colnames(sb))
#### delta_mu
boxplot(sb.plot$sb~sb.plot$delta.mu+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Average Bias",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(rmse)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(sb)-1),par("usr")[3]-0.03, labels = colnames(sb),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.mu),title=expression(delta[mu]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
#### delta_tau
boxplot(sb.plot$sb~sb.plot$delta.tau+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Average Bias",xlab="",xaxt="n",cex.axis=0.75,cex.lab=0.75)
axis(1, at=1.5 + 2*(1:ncol(sb)-1),labels=colnames(sb),tick=T,cex.axis=0.75,lwd.ticks=0)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.tau),title=expression(delta[tau]),
       fill=c("lightblue","tomato"),cex=15,ncol=2,bty="n")
#### tau
boxplot(sb.plot$sb~sb.plot$tau+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Average Bias",xlab="",xaxt="n",cex.axis=0.75,cex.lab=0.75)
axis(1, at=1.5 + 2*(1:ncol(sb)-1),labels=colnames(sb),tick=T,cex.axis=0.75,lwd.ticks=0)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$tau),title=expression(tau),
       fill=c("lightblue","tomato"),cex=15,ncol=2,bty="n")
#### Absolute bias
sb.plot <- reshape(cbind(params,abs(sb)),direction="long",varying=list(5:9),timevar="Model",times=colnames(sb),v.names = "sb")
sb.plot$Model <- factor(sb.plot$Model, levels=colnames(sb))
#### delta_mu
pdf("Figures/absolute_bias_delta_mu.pdf")
par(bty="L")
boxplot(sb.plot$sb~sb.plot$delta.mu+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Absolute Bias",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(sb)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(sb)-1),par("usr")[3]-0.03, labels = colnames(sb),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.mu),title=expression(delta[mu]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### delta_tau
pdf("Figures/absolute_bias_delta_tau.pdf")
par(bty="L")
boxplot(sb.plot$sb~sb.plot$delta.tau+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Absolute Bias",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(sb)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(sb)-1),par("usr")[3]-0.03, labels = colnames(sb),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.tau),title=expression(delta[tau]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### tau
pdf("Figures/absolute_bias_tau.pdf")
par(bty="L")
boxplot(sb.plot$sb~sb.plot$tau+sb.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Absolute Bias",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(sb)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(sb)-1),par("usr")[3]-0.03, labels = colnames(sb),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(sb)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$tau),title=expression(tau),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### Variance
v.plot <- reshape(cbind(params,v),direction="long",varying=list(5:9),timevar="Model",times=colnames(v),v.names = "v")
v.plot$Model <- factor(v.plot$Model, levels=colnames(v))
#### delta_mu
pdf("Figures/variance_delta_mu.pdf")
par(bty="L")
boxplot(v.plot$v~v.plot$delta.mu+v.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Variance",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(v)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(v)-1),par("usr")[3]-0.03, labels = colnames(v),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(v)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.mu),title=expression(delta[mu]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### delta_tau
pdf("Figures/variance_delta_tau.pdf")
par(bty="L")
boxplot(v.plot$v~v.plot$delta.tau+v.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Variance",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(v)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(v)-1),par("usr")[3]-0.03, labels = colnames(v),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(v)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$delta.tau),title=expression(delta[tau]),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
#### tau
pdf("Figures/variance_tau.pdf")
par(bty="L")
boxplot(v.plot$v~v.plot$tau+v.plot$Model,
        col=c("lightblue","tomato"),
        main="",ylab="Variance",xlab="",xaxt="n",cex.axis=1,cex.lab=1)
axis(1, at=1.5 + 2*(1:ncol(v)-1),labels=F,tick=T,lwd.ticks=0)
text(1.5 + 2*(1:ncol(v)-1),par("usr")[3]-0.03, labels = colnames(v),
     srt = 0, pos = 1, xpd = TRUE, cex=1)
for (i in 2.5 + 2*(1:(ncol(v)-1)-1)) abline(v=i,lty=2)
legend("topleft",legend=unique(params$tau),title=expression(tau),
       fill=c("lightblue","tomato"),cex=1,ncol=1,bty="n")
dev.off()
### CATE
boxplot(rmse.cate,ylab="RMSE (CATE)",cex.lab=0.75,cex.axis=0.75)
#### Bias
par(mfrow=c(2,2))
boxplot(sapply(bart.rdd.cate.bias,sapply,mean),main="BART-RDD",
        ylab="Bias",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
abline(h=0,lty=2)
####
# boxplot(sapply(bcf.cate.bias,sapply,mean),main="BCF",
#         ylab="Bias",cex.lab=0.75,cex.axis=0.75,xaxt="n",
#         xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
# axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
# abline(h=0,lty=2)
####
boxplot(sapply(sbart.cate.bias,sapply,mean),main="S-BART",
        ylab="Bias",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
abline(h=0,lty=2)
####
boxplot(sapply(tbart.cate.bias,sapply,mean),main="T-BART",
        ylab="Bias",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
abline(h=0,lty=2)
par(mfrow=c(1,1))
#### Variance
par(mfrow=c(2,2))
boxplot(do.call(cbind,bart.rdd.cate.var),main="BART-RDD",
        ylab="Variance",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
####
# boxplot(do.call(cbind,bcf.cate.var),main="BCF",
#         ylab="Variance",cex.lab=0.75,cex.axis=0.75,xaxt="n",
#         xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
# axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
####
boxplot(do.call(cbind,sbart.cate.var),main="S-BART",
        ylab="Variance",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
####
boxplot(do.call(cbind,tbart.cate.var),main="T-BART",
        ylab="Variance",cex.lab=0.75,cex.axis=0.75,xaxt="n",
        xlab=bquote("("~tau~","~delta[mu]~","~delta[tau]~")"))
axis(1, at=1:nrow(params),labels=ticks,cex.axis=0.55)
par(mfrow=c(1,1))
####
# test.obs <- 1
# cate.test <- data.frame(`BART-RDD`=rowMeans(bart$results[[test.obs]]),SBART=rowMeans(sbart$results[[test.obs]]),TBART=rowMeans(tbart$results[[test.obs]]),CATE=test$cate)
# par(mfrow=c(2,2))
# for (i in 1:4) plot(cate.test[,i],ylab=bquote(tau),main=colnames(cate.test)[i])
# par(mfrow=c(1,2))
# plot(cate.test$TBART~cate.test$BART.RDD,xlab="BART-RDD",ylab="T-BART",
#      cex.axis=0.75,cex.lab=0.75)
# abline(0,1,lty=2)
# plot(cate.test$CATE,cate.test[,1],pch=19,col=1,
#      xlim=c(min(cate.test$CATE),max(cate.test$CATE)+0.3),ylim=c(min(cate.test),max(cate.test)),
#      xlab="CATE",ylab="Models",xaxt="n",cex.axis=0.75,cex.lab=0.75)
# axis(1,at=cate.test$CATE+0.1,labels=round(cate.test$CATE,2),cex.axis=0.75)
# points(cate.test$CATE+0.1,cate.test[,2],pch=19,col=2)
# points(cate.test$CATE+0.2,cate.test[,3],pch=19,col=3)
# abline(0,1,lty=2)
# legend("topleft",legend=names(cate.test)[-4],col=1:3,pch=19,bty="n",cex=15,ncol=2)
# par(mfrow=c(1,1))
# ####
# pdf("Figures/sim_rmse.pdf")
# par(mfrow=c(2,2))
# matplot(rmse,type="n",bty="n",xlab=expression(kappa),ylab="RMSE",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(rmse[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,rmse[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,rmse[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(sqrt(sb),type="n",bty="n",xlab=expression(kappa),ylab="Absolute Bias",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(sqrt(sb)[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,sqrt(sb)[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,sqrt(sb)[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(sqrt(v),type="n",bty="n",xlab=expression(kappa),ylab="Std. Dev.",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(sqrt(v)[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,sqrt(v)[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,sqrt(v)[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# plot.new()
# legend("center",ncol = 2,col=1:6, lty=1, pch=20, cex=0.8, legend = names)
# dev.off()
# ###
# pdf("Figures/sim_intervals.pdf")
# par(mfrow=c(2,2))
# matplot(cov,type="n",bty="n",xlab=expression(kappa),ylab="Coverage",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(cov[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,cov[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,cov[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(ci,type="n",bty="n",xlab=expression(kappa),ylab="Interval size",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(ci[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,ci[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,ci[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(int,type="n",bty="n",xlab=expression(kappa),ylab="Distance to 95% coverage",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(int[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,int[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,int[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# plot.new()
# legend("center",ncol = 2,col=1:6, lty=1, pch=20, cex=0.8, legend = names)
# dev.off()
# ###
# pdf("Figures/sim_rmse_cate.pdf")
# par(mfrow=c(1,1))
# matplot(rmse.cate,type="n",bty="n",xlab=expression(kappa),ylab="RMSE",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(rmse.cate[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,rmse.cate[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,rmse.cate[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# legend("topleft",ncol = 2,col=1:4, lty=1, pch=20, cex=1, legend = names[1:4])
# dev.off()
# ###
# pdf("Figures/sim_intervals_cate.pdf")
# par(mfrow=c(2,2))
# matplot(rmse.cate,type="n",bty="n",xlab=expression(kappa),ylab="RMSE",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(rmse.cate[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,rmse.cate[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,rmse.cate[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(cov.cate,type="n",bty="n",xlab=expression(kappa),ylab="Coverage",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(cov.cate[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,cov.cate[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,cov.cate[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# matplot(ci.cate,type="n",bty="n",xlab=expression(kappa),ylab="Interval size",xaxt="n")
# axis(1,at=c(2,5,8),labels=c(0,1,2))
# matlines(ci.cate[1:3,],type="b",lty=1,pch=20)
# matlines(4:6,ci.cate[4:6,],type="b",lty=1,pch=20)
# matlines(7:9,ci.cate[7:9,],type="b",lty=1,pch=20)
# abline(v=3.5,lty=2)
# abline(v=6.5,lty=2)
# plot.new()
# legend("center",ncol = 2,col=1:4, lty=1, pch=20, cex=0.8, legend = names[1:4])
# dev.off()
# ###
# pdf("Figures/rmse_boxplot.pdf")
# par(mfrow=c(2,1))
# boxplot(rmse,main="ATE",bty="n",las=3,cex.axis=0.7)
# boxplot(rmse.cate,main="CATE",bty="n",las=3,cex.axis=0.7)
# dev.off()