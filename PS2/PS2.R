setwd("/Volumes/SSD Hans/Github/AdvMetrics/PS2")
library(tidyverse)

data <- read.table("data_ps2.out")
probit <- glm(V4~V2+V3, family = binomial(link=probit),data)
summary(probit)

th_bs <- read.table("theta_bs.txt")
se_bs <- read.table("std_bs.txt")
th_bs_u <- read.table("theta_bs_u.txt")
th_hat <- read.table("theta_hat.txt")

nom <- c("Alpha","Lambda","Gamma")
nom_se
names(th_bs) <- c("")


hist(th_bs[,1]-th_hat[1,1], breaks = 20, main = "")
hist(th_bs[,2]-th_hat[2,1], breaks = 20, main = "")
hist(th_bs[,3]-th_hat[3,1], breaks = 20, main = "")
# , main = "",1
hist(se_bs[,1]-th_hat[4,1], breaks = 20, main = "")
hist(se_bs[,2]-th_hat[5,1], breaks = 20, main = "")
hist(se_bs[,3]-th_hat[6,1], breaks = 20, main = "")

hist(th_bs[,1], breaks = 20, main = "")
hist(th_bs[,2], breaks = 20, main = "")
hist(th_bs[,3], breaks = 20, main = "")
hist(se_bs[,1], breaks = 20, main = "")
hist(se_bs[,2], breaks = 20, main = "")
hist(se_bs[,3], breaks = 20, main = "")

th_bs_h <- read.table("theta_bs.txt")
se_bs_h <- read.table("std_bs.txt")
th_bs_u_h <- read.table("theta_bs_u.txt")


plot(density(th_bs_h[,1]-th_hat[1,1]))

#Halton Estimates
# Theta (BS)=    6.39987465224453       0.121538915130835
# -3.91302110834853
# S.E. (BS)=   5.052107438816363E-002  2.376523950986178E-002

#Uniform estimates
# Theta =    6.38161772994278       0.120124370171031
# -3.90105580131420
# S.E. =   5.034238198792176E-002  2.376460453001593E-002  2.918537487836581E-002

results <- tibble( Variables = nom,
                   Coef = th_hat[1:3,1],
                   S.E  = th_hat[4:6,1],
                   `Coef (BS)` = th_bs_u[1:3,1],
                   `S.E. (BS)` = th_bs_u[4:6,1],
                   `Coef (BSH)` = th_bs_u_h[1:3,1],
                   `S.E. (BSH)` = th_bs_u_h[4:6,1])

res_m <- t(results)

save(results, res_m, probit, nom, file = "ps2.RData")
