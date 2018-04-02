##########################################
## Chetverikov, Kim, and Wilhelm (2018) ##
##           simulation setup           ##
##########################################

rm(list=ls())
library(foreach)
source("Test_SD_functions.R")

# basic setup
n     <- 100 # sample size
brep  <- 200
num.grid <- n

# Data generating process
X <- runif(n)
U <- rnorm(n, 0, 0.1)
Y1 <- U # under the null
Y2 <- -0.1*X + U   # under A1
Y3 <- -0.1*X^2 + U # under A2
Y4 <- 0.2*X - 0.2*exp(-250*(X-0.5)^2 ) + U # under A3
  
# save the results for each outcome
res1 <- CKW.test(Y1, X, brep, num.grid)
res2 <- CKW.test(Y2, X, brep, num.grid)
res3 <- CKW.test(Y3, X, brep, num.grid)
res4 <- CKW.test(Y4, X, brep, num.grid)

# display all the results on console
print( rbind(res1, res2, res3, res4) )
