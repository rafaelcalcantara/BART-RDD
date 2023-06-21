data <- readstata13::read.dta13("Downloads/sample_75_02.dta")
data <- data[data$endy %in% 1981:2001 & data$volquit!=1 & data$region!=0 & data$recall!=0,]
work <- readstata13::read.dta13("Downloads/work_history.dta")
data <- merge(data,work)

y <- log(data$ne_wage0)-log(data$wage0)
x <- data$duration
c = 3*365
w <- subset(data, select = c("age","female","bluecollar","austrian",
                             "education","married","iagrmining",
                             "icarsales","ihotel","imanufact",
                             "iservice","itransport","iwholesale",
                             "endmo","endy","region"))
sample <- complete.cases(cbind(y,w)) & is.finite(y)
y <- y[sample]
x <- x[sample]
w <- w[sample,]
saveRDS(data.frame(y=y,x=x,c=c,w=w),"~/Documents/Git/XBCF-RDD/Data/search.rds")
