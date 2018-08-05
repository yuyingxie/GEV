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
##########################
n1 = 200
p = 500

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
##Oracle 
Orc = solve(Sigma) %*% mu2

set.seed(6)
X1 =   mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
X2 =   mvrnorm(n = n1, mu = mu2, Sigma = Sigma)
X_tr = rbind(X1, X2)
Y_tr =c(rep(0, n1), rep(1, n1))
X1_test =   mvrnorm(n = n1 , mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n1 , mu = mu2, Sigma = Sigma)
X_test = rbind(X1_test, X2_test)
Y_test = c(rep(0,  n1), rep(1,  n1))
#set.seed(6)
mtotal = apply(X_tr, 2, mean)

cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas=c(1e-4,1e-3,1e-2,.1,1,10),lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda=cv.out$bestlambda, K=cv.out$bestK, lambda2=.3)
min( sum(out$ypred == (Y_test + 1) ), length(Y_test) - sum(out$ypred == (Y_test + 1) ) )

#res = FDACV(X_tr, Y_tr, fold = 3, lambda =seq(0.01, 0.1, length =30),  k = 1)
res = FDACV(X_tr, Y_tr, fold = 100, lambda = seq(0.02, 0.2, length =20),  k = 1)
lambda = res$lambdaopt[1]
#W3 = GEV_FISTA(Sigw, U_est, lambda = lambda, diff_thre = 1e-6, max_iter = 1000) # Same as FDA_GEV
res1 = rep(0, 20)
a = lambda = seq(0.05, 0.15, length =20)
for(i in 1:20){
    lambda = a[i]
    W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
    res1[i] = cor(Orc, W2)
}
cor(Orc, W2)
cor(Orc, obj$beta[-1])

Err.result = matrix(0, 50, 2)
for(i in 1:50){
    W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
    X1_test =   mvrnorm(n = n1 , mu = rep(0, p), Sigma = Sigma)
    X2_test =   mvrnorm(n = n1 , mu = mu2, Sigma = Sigma)
    X_test = rbind(X1_test, X2_test)
    Y_test = c(rep(0,  n1), rep(1,  n1))
    
    Ypredict =( (X_test - rep(1, dim(X_test)[1]) %*% t(mtotal)) %*% W2) > 0
    Err.result[i, 1] =min( sum(Ypredict == Y_test), sum(!Ypredict == Y_test)  )
    obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
    Err.result[i, 2] = obj$error * dim(X_test)[1]
  }



W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
cor(Orc, W2)
cor(Orc, obj$beta[-1])

X1_test =   mvrnorm(n = n1 , mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n1 , mu = mu2, Sigma = Sigma)
X_test = rbind(X1_test, X2_test)
Y_test = c(rep(0,  n1), rep(1,  n1))

result = matrix(0, 1, 3)
Ypredict =( (X_test - rep(1, dim(X_test)[1]) %*% t(mtotal)) %*% W2) > 0
result[1, 1] = min( sum(Ypredict == Y_test), sum(!Ypredict == Y_test)  )
obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
result[1, 2] = obj$error * dim(X_test)[1]
cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas=c(1e-4,1e-3,1e-2,.1,1,10),lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda=cv.out$bestlambda, K=cv.out$bestK, lambda2=.3)
result[1, 3] = min( sum(out$ypred == (Y_test + 1) ), length(Y_test) - sum(out$ypred == (Y_test + 1) ) )
result

plot( ( X_tr -  rep(1, length(Y_tr)) %*% t(mtotal))    %*% W2)
plot( cbind(1, X_tr) %*% obj$beta)
## DSDA 

