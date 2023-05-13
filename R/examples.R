setwd("~/Documents/Git/XBCF-RDD")
## Splitting constraints
x <- seq(-1,1,0.01)
y <- x
h <- 0.08
###
plot(x,y,"n")
abline(v=-h,lty=2)
abline(v=h,lty=2)
abline(v=0,lty=2)
segments(x0=-1,x1=-h-0.5,y0=-1,col="darkgreen")
segments(x0=-h-0.5,y0=-1,y1=-0.5,col="darkgreen")
segments(x0=-h-0.5,x1=-h-0.1,y0=-0.5,col="darkgreen")
segments(x0=-h-0.1,y0=-0.5,y1=0,col="darkgreen")
segments(x0=-h-0.1,x1=h+0.1,y0=0,col="darkgreen")
segments(x0=h+0.1,y0=0,y1=0.5,col="darkgreen")
segments(x0=h+0.1,x1=h+0.5,y0=0.5,col="darkgreen")
segments(x0=h+0.5,y0=0.5,y1=1,col="darkgreen")
segments(x0=h+0.5,x1=1,y0=1,col="darkgreen")
