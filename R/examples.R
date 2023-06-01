set.seed(0)
setwd("~/Documents/Git/XBCF-RDD")
## Example data
h <- 0.27
n <- 100
x <- runif(n,-1,1)
w <- runif(n,-1,1)
y <- w + x + (x>=0)*3 + rnorm(n)
###
temp <- function(tree,w,x,x0,y0,x1,y1)
{
    if (length(tree)==0) return()
    if (tree$var == 1)
    {
        segments(x0,w[tree$val],x1,w[tree$val])
        y0.l <- y0
        y1.l <- w[tree$val]
        y0.r <- w[tree$val]
        y1.r <- y1
        x0.l <- x0
        x0.r <- x0
        x1.l <- x1
        x1.r <- x1
    } else
    {
        segments(x[tree$val],y0,x[tree$val],y1)
        x0.l <- x0
        x1.l <- x[tree$val]
        x0.r <- x[tree$val]
        x1.r <- x1
        y0.l <- y0
        y0.r <- y0
        y1.l <- y1
        y1.r <- y1
    }
    temp(tree$left,w,x,x0.l,y0.l,x1.l,y1.l)
    temp(tree$right,w,x,x0.r,y0.r,x1.r,y1.r)
}
plot.trees <- function(w,x,tree,Owidth)
{
    plot(x,w,"p")
    temp(tree,w,x,par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4])
    abline(v=-Owidth,lty=2)
    abline(v=Owidth,lty=2)
    abline(v=0,lty=3)
}   
###
split <- function(tree,var,val)
{
    if (length(tree) != 0) stop("Not a leaf node")
    tree$var <- var
    tree$val <- val
    tree$left <- list()
    tree$right <- list()
    return(tree)
}
nbot <- function(tree)
{
    if (length(tree)==0) return(1)
    nbot(tree$left)+nbot(tree$right)
}
###
#### This tree is invalid because, although the partitions don't
#### cut the window, the middle partition has no points inside it
#### to the right of the cutoff. Besides, even if there were points
#### there, it's unlikely that this tree would be valid because of
#### the Opct condition. This condition ensures the partitions are
#### tighter around the window
root1 <- list()
root1 <- split(root1,1,order(w)[50])
root1$left <- split(root1$left,2,order(x)[25])
root1$right <- split(root1$right,1,order(w)[75])
root1$left$left <- split(root1$left$left,1,order(w)[20])
root1$left$right <- split(root1$left$right,2,order(x)[65])
root1$right$left <- split(root1$right$left,2,order(x)[20])
root1$right$right <- split(root1$right$right,2,order(x)[60])
###
root2 <- list()
root2 <- split(root2,1,order(w)[50])
root2$left <- split(root2$left,2,order(x)[25])
root2$right <- split(root2$right,1,order(w)[75])
root2$left$left <- split(root2$left$left,1,order(w)[20])
root2$left$right <- split(root2$left$right,2,order(x)[65])
root2$right$left <- split(root2$right$left,2,order(x)[20])
root2$right$right <- split(root2$right$right,2,order(x)[60])
###
root3 <- list()
root3 <- split(root3,1,order(w)[50])
root3$left <- split(root3$left,2,order(x)[20])
root3$right <- split(root3$right,2,order(x)[10])
root3$left$left <- split(root3$left$left,1,order(w)[20])
root3$left$right <- split(root3$left$right,2,order(x)[65])
root3$right$left <- split(root3$right$left,2,order(x)[25])
root3$right$right <- split(root3$right$right,2,order(x)[67])
###
par(mfrow=c(2,2))
plot.trees(w,x,root1,h)
plot.trees(w,x,root2,h)
plot.trees(w,x,root3,h)
plot.trees(w,x,root4,h)
par(mfrow=c(1,1))
