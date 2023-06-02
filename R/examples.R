set.seed(0)
setwd("~/Documents/Git/XBCF-RDD")
## Functions
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
    plot(x,w,"p",cex=0.75)
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
## Example data: w independent of x
h <- 0.25
n <- 100
x <- runif(n,-1,1)
w <- runif(n,-1,1)
y <- w + x + (x>=0)*3 + rnorm(n)
### Plots
#### This plot makes it clear that Owidth=0.25 is too narrow for
#### this dataset
plot(x,w)
abline(v=-h,lty=2)
abline(v=h,lty=2)
abline(v=0,lty=3)
#### This tree is invalid because it partitions the window
root1 <- list()
root1 <- split(root1,1,order(w)[50])
root1$left <- split(root1$left,2,order(x)[25])
root1$right <- split(root1$right,1,order(w)[75])
root1$left$left <- split(root1$left$left,1,order(w)[20])
root1$left$right <- split(root1$left$right,2,order(x)[65])
root1$right$left <- split(root1$right$left,2,order(x)[40])
root1$right$right <- split(root1$right$right,2,order(x)[60])
###
#### This tree is invalid because, although the partitions don't
#### cut the window, the middle partition has no points inside it
#### to the right of the cutoff. Besides, even if there were points
#### there, it's unlikely that this tree would be valid because of
#### the Opct condition. This condition ensures the partitions are
#### tighter around the window
root2 <- list()
root2 <- split(root2,1,order(w)[50])
root2$left <- split(root2$left,2,order(x)[25])
root2$right <- split(root2$right,1,order(w)[75])
root2$left$left <- split(root2$left$left,1,order(w)[20])
root2$left$right <- split(root2$left$right,2,order(x)[65])
root2$right$left <- split(root2$right$left,2,order(x)[20])
root2$right$right <- split(root2$right$right,2,order(x)[60])
###
#### This tree is invalid because, although the partitions don't
#### cut the window and the middle partition has points inside it
#### to the right of the cutoff, there are many more points
#### in the middle partition outside the window, violating the
#### Opct condition
root3 <- list()
root3 <- split(root3,1,order(w)[50])
root3$left <- split(root3$left,2,order(x)[20])
root3$right <- split(root3$right,2,order(x)[10])
root3$left$left <- split(root3$left$left,1,order(w)[20])
root3$left$right <- split(root3$left$right,2,order(x)[70])
root3$right$left <- split(root3$right$left,2,order(x)[25])
root3$right$right <- split(root3$right$right,2,order(x)[65])
###
root4 <- list()
root4 <- split(root4,1,order(w)[50])
root4$left <- split(root4$left,2,order(x)[20])
root4$right <- split(root4$right,2,order(x)[10])
root4$left$left <- split(root4$left$left,1,order(w)[20])
root4$left$right <- split(root4$left$right,2,order(x)[71])
root4$right$left <- split(root4$right$left,2,order(x)[24])
root4$right$right <- split(root4$right$right,2,order(x)[68])
###
png("Figures/bad_trees_1.png")
par(mfrow=c(1,3))
plot.trees(w,x,root1,h)
plot.trees(w,x,root2,h)
plot.trees(w,x,root3,h)
par(mfrow=c(1,1))
dev.off()
###
png("Figures/good_tree_1.png")
plot.trees(w,x,root4,0.36)
dev.off()
## Example data: w and x related
### In this setting, the partitions must reflect the relationship
### between w and x otherwise we end up with empty nodes.
### For example, in this case high w equals high x, meaning
### partitions at high values of w are less likely to feature
### x near the cutoff and thus inside the window. Similarly for
### low w and low x.
h <- 0.25
n <- 100
x <- runif(n,-1,1)
w <- sapply(x,function(i) rnorm(1,i,0.5))
y <- w + x + (x>=0)*3 + rnorm(n)
### Plots
#### This plot makes it clear that Owidth=0.25 is too narrow for
#### this dataset
plot(x,w)
abline(v=-h,lty=2)
abline(v=h,lty=2)
abline(v=0,lty=3)
#### This tree is invalid because it creates empty nodes in window
root1 <- list()
root1 <- split(root1,1,order(w)[49])
root1$left <- split(root1$left,2,order(x)[30])
root1$right <- split(root1$right,2,order(x)[32])
root1$left$left <- split(root1$left$left,1,order(w)[15])
root1$left$right <- split(root1$left$right,2,order(x)[65])
root1$right$left <- split(root1$right$left,2,order(x)[25])
root1$right$right <- split(root1$right$right,1,order(w)[85])
###
root2 <- list()
root2 <- split(root2,1,order(w)[49])
root2$left <- split(root2$left,2,order(x)[30])
root2$right <- split(root2$right,2,order(x)[32])
root2$left$left <- split(root2$left$left,1,order(w)[15])
root2$left$right <- split(root2$left$right,2,order(x)[65])
root2$right$left <- split(root2$right$left,2,order(x)[25])
root2$right$right <- split(root2$right$right,2,order(x)[67])
###
png("Figures/trees_2.png")
par(mfrow=c(1,2))
plot.trees(w,x,root1,h)
plot.trees(w,x,root2,h)
par(mfrow=c(1,1))
dev.off()
