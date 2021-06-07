#### Setting ####

setwd("/Volumes/SSD Hans/Github/AdvMetrics/Final")
library(tidyverse)
library(sandwich)
library(ggplot2)
data <- read_csv("res.csv")

#### Q3-Q4 testing routines graphs ----
par(mfrow=c(3,1))
plot(NA, xlim=range(data$A2), ylim=range(c(data$V2_0, data$V2_1)),
     xlab=expression(A[2]), ylab=expression(V[2]))
lines(data$V2_0~data$A2, lty=1, lwd=4, col=coul[1])
lines(data$V2_1~data$A2, lty=1, lwd=4, col=coul[2])
title("t=2")
legend("bottomright", legend=c(expression(L[1]==0), expression(L[1]==1)),
       lty = 1, col = coul, lwd = 2)
abline(v=1.1, col = coul[1], lty=2)
abline(v=1.37, col = coul[2], lty=2)

plot(NA, xlim=range(data$A2), ylim=range(c(data$A3_0, data$A3_1)),
     xlab=expression(A[2]), ylab=expression(A[3]))
lines(data$A3_0~data$A2, lty=1, lwd=4, col=coul[1])
lines(data$A3_1~data$A2, lty=1, lwd=4, col=coul[2])
title("t=2")
legend("bottomright", legend=c(expression(L[1]==0), expression(L[1]==1)),
       lty = 1, col = coul, lwd = 2)
abline(v=1.1, col = coul[1], lty=2)
abline(v=1.37, col = coul[2], lty=2)

plot(NA, xlim=range(data$A2), ylim=range(c(data$L2_0, data$L2_1)),
     xlab=expression(A[2]), ylab=expression(L[2]))
lines(data$L2_0~data$A2, lty=1, lwd=4, col=coul[1])
lines(data$L2_1~data$A2, lty=1, lwd=4, col=coul[2])
title("t=2")
legend("bottomright", legend=c(expression(L[1]==0), expression(L[1]==1)),
       lty = 1, col = coul, lwd = 2)
abline(v=1.1, col = coul[1], lty=2)
abline(v=1.37, col = coul[2], lty=2)

par(mfrow=c(3,1))
plot(NA, xlim=range(res_t1$A1), ylim=range(res_t1$V1),
     xlab=expression(A[1]), ylab=expression(V[1]))
lines(res_t1$V1~res_t1$A1, lty=1, lwd=4, col=coul)
title("t=1")
# legend("bottomright", legend=c(expression(L[1]==0), expression(L[1]==1)),
#        lty = 1, col = c("red",coul), lwd = 2)
abline(v=1.71, col = coul, lty=2)

plot(NA, xlim=range(res_t1$A1), ylim=range(res_t1$A2),
     xlab=expression(A[1]), ylab=expression(A[2]))
lines(res_t1$A2~res_t1$A1, lty=1, lwd=4, col=coul)
title("t=1")
abline(v=1.71, col = coul, lty=2)

plot(NA, xlim=range(res_t1$A1), ylim=range(res_t1$L1),
     xlab=expression(A[1]), ylab=expression(L[1]))
lines(res_t1$L1~res_t1$A1, lty=1, lwd=4, col=coul)
title("t=1")
abline(v=1.71, col = coul, lty=2)


### Check Fortran regs ----
df <- read.table("data.out")
names(df)=c("A1","L1","Y1","A2","L2","Y2","A3")


summary(lm(A2~Y1+A1+L1, df[1:1000,]))
var(residuals(lm(A2~Y1+A1+L1, df[1:1000,])))

summary(lm(A3~Y2+A2+L2, df[1:1000,]))
var(residuals(lm(A3~Y2+A2+L2, df[1:1000,])))

df2=df[1:1000,]

pick=df2$L1==1
mean(log(df2[pick,"Y1"]))
var(log(df2[pick,"Y1"]))


pick2=df2$L2==1
summary(lm(log(Y2)~L1, df2[pick2,]))
var(residuals(lm(log(Y2)~L1, df2[pick2,])))

sum(pick2)

#### Indirect inference ----

iid <- read_csv("fake_df_mpi.csv", skip = 9, col_names = F)
names(iid)=c("A1","L1","Y1","A2","L2","Y2","A3")
coul <- c("#4f2683","#807f83","#14bdeb","#00100b","#a6d49f")

## Distribution of assets

### Base R
hist(iid$A1)
plot(iid$A1,iid$A2)
hist(iid$A2)


par(mfrow=c(1,1))

plot(NA, xlim=c(0,5),
     ylim=range(c(density(iid$A1)$y, density(iid$A2)$y, density(iid$A3)$y)),
     xlab="Assets", ylab="")
lines(density(iid$A1), lwd = 2, col = coul[1])
lines(density(iid$A2), lwd = 2, col = coul[2])
lines(density(iid$A3), lwd = 2, col = coul[3])
title("Distribution of asssets")
legend("topright", legend = c(expression(A[1]), expression(A[2]), expression(A[3])),
       lty = 1, col = coul, lwd = 2)

### GGPLOT

dens <- iid %>%
        pivot_longer(starts_with("A"),
                     names_to = "Period",
                     values_to = "assets" )

dens %>% ggplot(aes(x=assets, group=period, fill=period))+
        geom_density(alpha=0.4, adjust=1.5)+
        xlim(0,5)+
        theme_classic()+
        theme(legend.title=element_blank())+
        scale_fill_manual(values = coul)+
        xlab("") + ylab("")+ggtitle("Asset distribution")

## unemployment rate
## example base r
emp <- data.frame( t1=c((sum(iid$L1)/1008)*100, (1-sum(iid$L1)/1008)*100),
                       t2=c((sum(iid$L2)/1008)*100, (1-sum(iid$L2)/1008)*100))
rownames(emp) <- c("Emp","Unemp")

barplot(as.matrix(emp), col=coul, border="white", main = "Employment rate")

## Example ggplot
iid %>% pivot_longer(starts_with("L"),
                     names_to = "Period",
                     values_to = "Labour") %>%
        mutate(value=1, Labour=as.factor(Labour)) %>%
        arrange(desc(Labour)) %>%
        ggplot(aes(x=Period, fill=Labour, y=value))+
        geom_bar(position = "stack", stat = "identity")+
        theme_classic()+
        scale_fill_manual(values = coul[c(2,1)])+
        theme(legend.position = "none")+
        xlab("") + ylab("")+ggtitle("Employment rate")


## table of resutls
iivt <- read_csv("ev.csv")


### counterfactuals ----
imdf <- read_csv("fake_df_mpi_6a.csv", skip = 9, col_names = F)
names(imdf)=c("A1","L1","Y1","A2","L2","Y2","A3")
irdf <- read_csv("fake_df_mpi_6b.csv", skip = 9, col_names = F)
names(irdf)=c("A1","L1","Y1","A2","L2","Y2","A3")

# Want to compare three counterfactual exercises in one graph

empdf <- iid %>% pivot_longer(starts_with("L"),
                     names_to = "Period",
                     values_to = "Labour") %>%
        mutate(value=1, Labour=as.factor(Labour), CF="Base case") %>%
        select(Period, Labour, value, CF) %>%
        arrange(desc(Labour))

tdf1 <- imdf %>% pivot_longer(starts_with("L"),
                              names_to = "Period",
                              values_to = "Labour") %>%
        mutate(value=1, Labour=as.factor(Labour), CF="I min increase") %>%
        select(Period, Labour, value, CF) %>%
        arrange(desc(Labour))

tdf2 <- irdf %>% pivot_longer(starts_with("L"),
                              names_to = "Period",
                              values_to = "Labour") %>%
        mutate(value=1, Labour=as.factor(Labour), CF="r increase") %>%
        select(Period, Labour, value, CF) %>%
        arrange(desc(Labour))

empdf <- rbind(empdf, tdf1, tdf2)
empdf %>% filter(Period=="L1")  %>%
        mutate(CF=as.factor(CF)) %>%
        ggplot(aes(x=CF, fill=Labour, y=value))+
        geom_bar(position = "stack", stat = "identity")+
        theme_classic()+
        scale_fill_manual(values = coul[c(2,1,3,4,5,6)])+
        theme(legend.position = "none")+
        xlab("") + ylab("")+ggtitle("Employment rate t=1")

empdf %>% filter(Period=="L2")  %>%
        mutate(CF=as.factor(CF)) %>%
        ggplot(aes(x=CF, fill=Labour, y=value))+
        geom_bar(position = "stack", stat = "identity")+
        theme_classic()+
        scale_fill_manual(values = coul[c(2,1,3,4,5,6)])+
        theme(legend.position = "none")+
        xlab("") + ylab("")+ggtitle("Employment rate t=2")

### Ass distribution

asdf <- iid %>% pivot_longer(starts_with("A"),
                     names_to = "Assets",
                     values_to = "Values" ) %>%
        mutate(CF="Base case") %>%
        select(Assets, Values, CF)

dft1<- imdf %>% pivot_longer(starts_with("A"),
                             names_to = "Assets",
                             values_to = "Values" ) %>%
        mutate(CF="I min increase") %>%
        select(Assets, Values, CF)

dft2<- irdf %>% pivot_longer(starts_with("A"),
                             names_to = "Assets",
                             values_to = "Values" ) %>%
        mutate(CF="r increase") %>%
        select(Assets, Values, CF)

asdf <- rbind(asdf,dft1, dft2)


asdf %>% filter(Assets=="A3") %>%
        mutate(CF=as.factor(CF)) %>%
        ggplot(aes(x=Values, group=CF, fill=CF, color=CF))+
        geom_density(alpha=0.2)+
        xlim(0,3)+
        theme_classic()+
        theme(legend.title=element_blank())+
        scale_fill_manual(values = coul)+
        scale_color_manual(values = coul)+
        xlab("") + ylab("")+ggtitle("Asset distribution A3")

asdf %>% filter(Assets=="A2") %>%
        mutate(CF=as.factor(CF)) %>%
        ggplot(aes(x=Values, group=CF, fill=CF, color=CF))+
        geom_density(alpha=0.2)+
        xlim(0,3)+
        theme_classic()+
        theme(legend.title=element_blank())+
        scale_fill_manual(values = coul)+
        scale_color_manual(values = coul)+
        xlab("") + ylab("")+ggtitle("Asset distribution A2")


asdf

### saving results to read in rmarkdown ----
save(iivt, file = "final.RData")