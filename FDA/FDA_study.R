library(PMA)
library(MASS)
library(RSpectra)
library(msda)
library(dsda)
library(penalizedLDA)
source("../GEV.R")

args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    case.id = 1  
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

if(Type == "M"){
    set.seed(case.id)
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
    set.seed(case.id)
    X1 =   mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
    X2 =   mvrnorm(n = n, mu = mu2, Sigma = Sigma)
    X3 =   mvrnorm(n = n, mu = mu3, Sigma = Sigma)
    
    X_tr = rbind(X1, X2, X3)
    Y_tr =c(rep(0, n), rep(1, n), rep(2, n))
    mtotal = apply(X_tr, 2, mean)
    
    n2 = 2 * n
    X1_test =   mvrnorm(n = n2, mu = rep(0, p), Sigma = Sigma)
    X2_test =   mvrnorm(n = n2 , mu = mu2, Sigma = Sigma)
    X3_test =   mvrnorm(n = n2 , mu = mu3, Sigma = Sigma)
    X_test = rbind(X1_test, X2_test, X3_test)
    Y_test = c(rep(0,  n2 ), rep(1,  n2), rep(2,  n2))
}else{    
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
    set.seed(case.id)
    X1 =   mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma) # n x p: 200 x 50
    X2 =   mvrnorm(n = n, mu = mu2, Sigma = Sigma)
    X_tr = rbind(X1, X2)
    Y_tr =c(rep(0, n), rep(1, n))
    X1_test =   mvrnorm(n = n * 2, mu = rep(0, p), Sigma = Sigma)
    X2_test =   mvrnorm(n = n * 2 , mu = mu2, Sigma = Sigma)
    X_test = rbind(X1_test, X2_test)
    Y_test = c(rep(0,  n * 2), rep(1,  n * 2)) 
}

result = matrix(0, 1, 4)

res1 = FDACV(X_tr, Y_tr, fold = 6, lambda = seq(0.1, 0.4, length = 20),  k = 2, max_iter = 2000)
lambda = res1$lambdaopt[1]

FDA_result = FDA_pred(X_tr, Y_tr, X_test, Y_test, lambda = lambda,
        diff_thre = 1e-6, max_iter = 3000,  standardize = TRUE, k = 2 )
result[1, 1] = FDA_result$error

if(Type == "M"){
    Msda.cv = cv.msda(x = X_tr, y = Y_tr + 1, nfolds = 5, lambda.opt = "max")
    lambda.min = Msda.cv$lambda.min
    id.min = which(Msda.cv$lambda == lambda.min)
    Msda_pred = predict(Msda.cv$msda.fit, X_test)[ , id.min]
    result[1, 2] = sum(Msda_pred != (Y_test + 1))
}else{    
    obj<-dsda(X_tr, Y_tr, X_test, Y_test)  ##perform direct sparse discriminant analysis
    result[1, 2] = obj$error * dim(X_test)[1]
}

cv.out =  PenalizedLDA.cv(X_tr, Y_tr + 1, lambdas = c(1e-4, 1e-3, 1e-2, .1, 1, 10), lambda2=.3)
out = PenalizedLDA(X_tr, Y_tr + 1, xte = X_test,  lambda = cv.out$bestlambda, K = cv.out$bestK, lambda2 = .3)
if(Type == "M"){
    result[1, 3] = sum(out$ypred[, 2] != (Y_test + 1))
}else{
    result[1, 3] = sum(out$ypred != (Y_test + 1))
}
result[1, 4] = Orc_pred(mu, Sigma, X_test, Y_test)$error
res = list(result = result, FDA_result = FDA_result, cv_res = res1)
save(res,  file = paste("./Result/", Type, "_p", p, "_n_", n,"_id", case.id,"_Res.Rdata", sep = "")) 
