# TODO: Add comment
# 
# Author: xyy
###############################################################################
library(PMA)
library(MASS)
library(RSpectra)
library(msda)
library(dsda)
library(penalizedLDA)
source("../GEV.R")

##########################
#### Binary Case ############
n1 = 250
p = 1000 

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
#Omega = NearestNeighborSigma1(p, 5, 10)$Omega
#Sigma = solve(Omega)

set.seed(6)
X1 =   mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
X2 =   mvrnorm(n = n1, mu = mu2, Sigma = Sigma)
X_tr = rbind(X1, X2)
Y_tr =c(rep(0, n1), rep(1, n1))
#set.seed(6)
mtotal = apply(X_tr, 2, mean)

res = FDACV(X_tr, Y_tr, fold = 3, lambda =seq(0.001, 0.1, length =30),  k = 1)

lambda = res$lambdaopt[1]
#W2 = GEV_FISTA(Sigw, U_est, lambda = lambda, diff_thre = 1e-6, max_iter = 1000)
W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE )

X1_test =   mvrnorm(n = n1 , mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n1 , mu = mu2, Sigma = Sigma)
X_test = rbind(X1_test, X2_test)
Y_test = c(rep(0,  n1), rep(1,  n1))

obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
obj$error * dim(X_test)[1]
Ypredict =( (X_test - rep(1, dim(X_test)[1]) %*% t(mtotal)) %*% W2) > 0
min( sum(Ypredict == Y_test), sum(!Ypredict == Y_test)  )

plot( ( X_tr -  rep(1, n) %*% t(mtotal))    %*% W2, xlim = c(0, 400))
plot( cbind(1, X_tr) %*% obj$beta, xlim = c(0, 400))

## DSDA 



##Oracle 
p = 10
n = 100
X =    mvrnorm(n = n, mu = rep(0,  p), Sigma = diag(p))
Y = X %*% rep(1, p) + rnorm(n, mean = 0, sd = 0.1)
res = lm(Y~ X)
sum(res$residuals^2)

p = 10
n = 100
x = rep(0, p)
A = diag(p)
Y = matrix(0, n, p)
for(i in 1:n){
    Y[i, ] = A %*% x + rnorm(p, 0, 0.1) 
}
