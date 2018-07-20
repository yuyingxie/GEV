# TODO: Add comment
# 
# Author: xyy
###############################################################################


library(PMA)
library(MASS)
library(RSpectra)
library(msda)
source("../GEV.R")

##########################
#### Binary Case ############
n = 200
p = 50 

mu2 = rep(0, p)
mu2[(1:20) * 2] = 0.5
Sigma = diag(p)

for(m in 1:5){
    mm = p/5 * (m - 1)
    for(i in 1:(p/5)){
        for(j in 1:(p/5)){
            Sigma[i + mm, j + mm] = 0.7^{abs(i - j)}
        }
    }
}

X1 =   mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
X2 =   mvrnorm(n = n, mu = mu2, Sigma = Sigma)

m1 = apply(X1, 2, mean)
m2 = apply(X2, 2, mean)
X1c = X1 - rep(1, n) %*% t(m1) 
X2c = X2 - rep(1, n) %*% t(m2) 
Sigw = (t(X1c ) %*% X1c +  t(X2c) %*% X2c) / (2 * n)
Sigb_est = (n * m1 %*% t(m1) + n * m2 %*% t(m2))/ (2 * n )  
Sigb = (n * 0 + n * mu2 %*% t(mu2))/ (2 * n )  

X1_test =   mvrnorm(n = n * 2, mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n * 2, mu = mu2, Sigma = Sigma)

X_tr = rbind(X1, X2)
Y_tr =c(rep(1, n), rep(2, n))

obj.cv = cv.msda(X_tr, Y_tr)
## MSDA
        