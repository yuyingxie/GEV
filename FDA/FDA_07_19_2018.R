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
n1 = 150
p = 500
set.seed(6)
mu2 = rep(0, p)
mu2[(1:20) ] = rnorm(20, 0.3, 0.5)
Sigma = diag(p)
mu = cbind(rep(0, p), mu2)

for(m in 1:5){
    mm = p/5 * (m - 1)
    for(i in 1:(p/5)){
        for(j in 1:(p/5)){
            Sigma[i + mm, j + mm] = 0.7^{abs(i - j)}
        }
    }
}
#Omega = NearestNeighborSigma1(p, 5, 6)$Omega
#Sigma = solve(Omega)
##Oracle 
Orc = solve(Sigma) %*% mu2

set.seed(36)
X1 =   mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
X2 =   mvrnorm(n = n1, mu = mu2, Sigma = Sigma)
X_tr = rbind(X1, X2)
Y_tr =c(rep(0, n1), rep(1, n1))
X1_test =   mvrnorm(n = n1 * 2, mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n1 * 2, mu = mu2, Sigma = Sigma)
X_test = rbind(X1_test, X2_test)
Y_test = c(rep(0,  n1 * 2), rep(1,  n1 * 2))
#set.seed(6)

result = matrix(0, 1, 4)
#res = FDACV(X_tr, Y_tr, fold = 3, lambda =seq(0.01, 0.1, length =30),  k = 1)
res = FDACV(X_tr, Y_tr, fold = 6, lambda = seq(0.1, 0.3, length =15),  k = 1)
lambda = res$lambdaopt[1]
FDA_result = FDA_pred(X_tr, Y_tr, X_test, Y_test, lambda = res$lambdaopt[1],
        diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 2 )
result[1, 1] = FDA_result$error

obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
result[1, 2] = obj$error * dim(X_test)[1]

cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas = c(1e-4, 1e-3, 1e-2, .1, 1, 10), lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda = cv.out$bestlambda, K = cv.out$bestK, lambda2 = .3)
result[1, 3] = sum(out$ypred != (Y_test + 1))

result[1, 4] = Orc_pred(mu, Sigma, X_test, Y_test)$error
result


#########################
#########################
Res = matrix(0, 50, 4)

type = "M"
p = 500
n = 150

for(i in 1:50){
    load(paste(type, "_p", p, "_n_", n, "_id", i, "_Res.Rdata", sep = ''))
    Res[i, ] = res$result
}

apply(Res, 2, mean)


result = matrix(0, 50, 4)
for(i in 1:50){
    load("/Users/yuyingxie/Dropbox (ValdarLab)/Michigan_State_Idea/Manuscript/GEV/GEV/FDA/Result/B_p500_n_150_id1_Res.Rdata")
}
















obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
cor(Orc, obj$beta[-1])
cor(Orc, W2)
#W3 = GEV_FISTA(Sigw, U_est, lambda = lambda, diff_thre = 1e-6, max_iter = 1000) # Same as FDA_GEV
res1 = rep(0, 20)
a = lambda = seq(0.05, 0.10, length =20)
for(i in 1:20){
    lambda = a[i]
    W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
    res1[i] = cor(Orc, W2)
}
obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
cor(Orc, obj$beta[-1])
cor(Orc, W2)

Err.result = matrix(0, 50, 3)
for(i in 1:50){
    n2 = n1 * 2
    W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
    X1_test =   mvrnorm(n = n2 , mu = rep(0, p), Sigma = Sigma)
    X2_test =   mvrnorm(n = n2 , mu = mu2, Sigma = Sigma)
    X_test = rbind(X1_test, X2_test)
    Y_test = c(rep(0,  n2), rep(1,  n2))
    
    Ypredict =( (X_test - rep(1, dim(X_test)[1]) %*% t(mtotal)) %*% W2) > 0
    Err.result[i, 1] = min( sum(Ypredict == Y_test), sum(!Ypredict == Y_test)  )
    obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
    Err.result[i, 2] = obj$error * dim(X_test)[1]
    cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas=c(1e-4,1e-3,1e-2,.1,1,10),lambda2=.3)
    out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda=cv.out$bestlambda, K=cv.out$bestK, lambda2=.3)
    Err.result[i, 3] = min( sum(out$ypred == (Y_test + 1) ), length(Y_test) - sum(out$ypred == (Y_test + 1) ) )
  }

W2 = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 1 )
cor(Orc, W2)
cor(Orc, obj$beta[-1])
n2 = 2 * n1
X1_test =   mvrnorm(n = n2, mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n2 , mu = mu2, Sigma = Sigma)
X_test = rbind(X1_test, X2_test)
Y_test = c(rep(0,  n2 ), rep(1,  n2))

result = matrix(0, 1, 3)
Ypredict =( (X_test - rep(1, dim(X_test)[1]) %*% t(mtotal)) %*% W2) > 0
result[1, 1] = min( sum(Ypredict == Y_test), sum(!Ypredict == Y_test)  )
obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
result[1, 2] = obj$error * dim(X_test)[1]
cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas=c(1e-4,1e-3,1e-2,.1,1,10),lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda=cv.out$bestlambda, K=cv.out$bestK, lambda2=.3)
result[1, 3] = min( sum(out$ypred == (Y_test + 1) ), length(Y_test) - sum(out$ypred == (Y_test + 1) ) )
result


FDA_pred(X_tr, Y_tr, X_test, Y_test,lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 2 )
plot( ( X_tr -  rep(1, length(Y_tr)) %*% t(mtotal))    %*% W2)
plot( cbind(1, X_tr) %*% obj$beta)
## DSDA 


##############
## K = 3
###############
n1 = 150
p = 500
set.seed(6)
mu2 = mu3 = rep(0, p)
mu2[(1:20) ] = rnorm(20, 0.3, 0.5)
mu3[(21:40) ] = -rnorm(20, -0.5, 0.5)
mu = cbind(rep(0, p), mu2, mu3)

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
#Orc = solve(Sigma) %*% mu2

set.seed(6) # 36 is good
X1 =   mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
X2 =   mvrnorm(n = n1, mu = mu2, Sigma = Sigma)
X3 =   mvrnorm(n = n1, mu = mu3, Sigma = Sigma)

X_tr = rbind(X1, X2, X3)
Y_tr =c(rep(0, n1), rep(1, n1), rep(2, n1))
mtotal = apply(X_tr, 2, mean)

n2 = 2 * n1
set.seed(12345)
X1_test =   mvrnorm(n = n2, mu = rep(0, p), Sigma = Sigma)
X2_test =   mvrnorm(n = n2 , mu = mu2, Sigma = Sigma)
X3_test =   mvrnorm(n = n2 , mu = mu3, Sigma = Sigma)
X_test = rbind(X1_test, X2_test, X3_test)
Y_test = c(rep(0,  n2 ), rep(1,  n2), rep(2,  n2))

result = matrix(0, 1, 4)

res = FDACV(X_tr, Y_tr, fold = 6, lambda = seq(0.1, 0.3, length = 10),  k = 2, max_iter = 2000)
lambda = res$lambdaopt[1]

FDA_result = FDA_pred(X_tr, Y_tr, X_test, Y_test, lambda = lambda,
        diff_thre = 1e-6, max_iter = 3000,  standardize = TRUE, k = 2 )
result[1, 1] = FDA_result$error

Msda.cv = cv.msda(x = X_tr, y = Y_tr + 1, nfolds = 5, lambda.opt = "max")
lambda.min = Msda.cv$lambda.min
id.min = which(Msda.cv$lambda == lambda.min)
Msda_pred = predict(Msda.cv$msda.fit, X_test)[ , id.min]
result[1, 2] = sum(Msda_pred != (Y_test + 1))

cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas = c(1e-4, 1e-3, 1e-2, .1, 1, 10), lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda = cv.out$bestlambda, K = cv.out$bestK, lambda2 = .3)
result[1, 3] = sum(out$ypred[, 2] != (Y_test + 1))

result[1, 4] = Orc_pred(mu, Sigma, X_test, Y_test)$error
result

a = rep(0, 10)
lam = seq(0.02, 0.1, length = 10)
for(i in 1:10){
    FDA_result = FDA_pred(X_tr, Y_tr, X_test, Y_test, lambda = lam[i],
            diff_thre = 1e-6, max_iter = 3000,  standardize = TRUE, k = 2 )
    a[i] = FDA_result$error
}





Res = matrix(0, 50, 4)
type = "M"
p = 500
n = 100

for(i in 1:50){
    load(paste(type, "_p", p, "_n_", n, "_id", i, "_Res.Rdata", sep = ''))
    Res[i, ] = res$result
}

apply(Res, 2, mean)





W = FDA_GEV(X_tr, Y_tr, lambda = lambda, diff_thre = 1e-6, max_iter = 1000,  standardize = TRUE, k = 2)
W_sc = 1/ sqrt(diag(t(W) %*% W))
W = W %*% diag(W_sc)
Orc = solve(Sigma) %*% mu2
Orc = cbind(Orc, solve(Sigma) %*% mu3)

plot( ( X_tr -  rep(1, length(Y_tr)) %*% t(mtotal))    %*% W)
plot( ( X_tr -  rep(1, length(Y_tr)) %*% t(mtotal))    %*% Orc, xlim = c(-20, 30), ylim = c(-20, 30))
