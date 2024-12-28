polynomial <- function(x,a,b,c,d) a*x^3+b*x^2+c*x+d
x <- seq(-0.1,0.1,length.out=1000)
plot(x,polynomial(x,10,1,0,1))