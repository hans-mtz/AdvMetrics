setwd("/Volumes/SSD Hans/Github/AdvMetrics/PS3")
library(tidyverse)
library(sandwich)
data <- read.table("data_ps2.out")

ols <- lm(V4~V2+V3, data)

vcovHC(ols,type = "HC0")
vcov(ols)
summary(ols)


x <- cbind(1,data[,2],data[,3])
x
xx <- t(x)%*%x
xx
solve(xx)
e <- residuals(ols)
e_s <- sum(e^2)/50000
# xux <- 0
# for (i in 1:50000) {
#   xux[i] <- (x[i,]%*%t(x[i,]))*e_s[i]
# }
# xux <- sum(xux)
varm<- 50000*solve(xx)%*%(xx*e_s)%*%solve(xx)
varm

solve(xx)*e_s
