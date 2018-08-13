# TODO: Add comment
# 
# Author: xyy
###############################################################################

library(psych)

Get_U = function(Omega, K){
    # Get U out of Omega
    # Input:
    # 			Omega is a p x p matrix
    #          K is the selected rank
    #temp.eigen = eigen(Omega)
    temp = eigs_sym(Omega, K * 2)
    #U = temp.eigen$vector[, 1:K] %*% diag(sqrt(temp.eigen$value[1:K]) )    
    if(K == 1){
        U = as.matrix(temp$vector[, 1] * sqrt(temp$value[1:K]) )
    }else{
    U = temp$vector[, 1:K] %*% diag(sqrt(temp$value[1:K]) )    
    }
    return(U)
}

Get_U_CCA = function(Omega, K){
    # Get U out of Omega
    # Input:
    # 			Omega is a p x p matrix
    #          K is the selected rank for CCA, we will pick 2K of the eigen vector of Ome since it block with zero
    #			in the diag	            
    #temp.eigen = eigen(Omega)
    temp = eigs_sym(Omega, K * 2)
    #U = temp.eigen$vector[, 1:K] %*% diag(sqrt(temp.eigen$value[1:K]) )    
    U = temp$vector[, 1:(2*K)] %*% diag(sqrt(temp$value[1:(2 *K)]) )    
    return(U)
}

GEV_FISTA = function(Sigma, U, lambda = 0.1, diff_thre = 2e-5, max_iter = 1000){
    #####################################################
    ##### Sigma is a p x p matrix  can be Cov(X) or Cov(X - E(X|Y))   #####
    ##### U is the p x k  matrix s.t.   Cov(E(X| Y)) = U U^T                    #####
    ##### lambda is the tuning parameter  controlling sparsity            #####
    ##### Objective function: tr( t(W) %*% \Sigma %*% W) * 0.5 - tr(t(W) %*%U) + lambda *  sum (abs(W))
    ######################################################
    k = dim(U)[2]
    d = dim(U)[1]
    #W = matrix(0.1, d, k)
    W = diag(1/ diag(Sigma) ) %*% U
    dif_value = 1
    count = 1
    W_new = W_y_new = W
    W_y = W
    eig1 = eigen(Sigma)$value[1]
    alpha = 1/eig1   # step size
    t.0 = 1 
    temp = alpha * (Sigma %*% W  - U)
    
    while((dif_value > diff_thre) & (count < max_iter)){
        t_new = 0.5 * (1 + sqrt(1 + 4 * t.0^2))
        for(j in 1:k){
            # D.L = D.Huber(Sigma[j, ] %*% W[, j] - U[, j], a = a)
            # W_new[, j] = sapply(W_y[, j] - temp[, j],  Shrink, lambda = lambda * alpha)
            W_new[, j] = mapply(Shrink, W_y[, j] - temp[, j],  lambda * alpha)
            W_y[, j] = W_new[, j] + (t.0 - 1)/ t_new * (W_new[, j] - W[, j])            # W_y is W_X in the manu
                                                                   # W_new is W_Y^{(k + 1)} in the manu and W is W_Y^{(k)}
        }
      
        temp = alpha * (Sigma %*% W_y  - U)
        dif_value = sum(abs(W - W_new))/sum(abs(W))
        W = W_new
        t.0 = t_new
        count = count + 1
        
        res = tr( t(W) %*% Sigma %*% W) * 0.5 - tr(t(W) %*%U) + sum(lambda * apply(abs(W), 1, sum)) 
        if(sum(abs(W)) == 0){
            break
        }
        cat(sep = "", "[",count,"]","dif=", round(dif_value, 5), "; Loss=", res, "\n")
        
    }
    return(W)
}

GEV_ISTA = function(Sigma, U, lambda = 0.1, diff_thre = 2e-5, max_iter = 1000){
    #####################################################
    ##### X is a n x p matrix  #####
    ##### Y is a n x m             #####
    ##### lambda is the tuning parameter  controlling sparsity            #####
    #############################################
    d = dim(U)[2]
    p = dim(U)[1]
   # W = matrix(0.1, p, d)
    W = diag(1/ diag(Sigma) ) %*% U
    dif_value = 1
    count = 1
    W_new = W_y_new = W
    W_y = W
    eig1 = eigen(Sigma)$value[1]
    alpha = 1/eig1   # step size
    t.0 = 1 
    temp = alpha * (Sigma %*% W  - U)
    while((dif_value > diff_thre) & (count < max_iter)){
        t_new = 0.5 * (1 + sqrt(1 + 4 * t.0^2))
        for(j in 1:d){
            # D.L = D.Huber(Sigma[j, ] %*% W[, j] - U[, j], a = a)
            W_new[, j] = sapply(W_y[, j] - temp[, j],  Shrink, lambda = lambda * alpha)
            W_y[, j] = W_new[, j]         
        }
        
        temp = alpha * (Sigma %*% W_y  - U)
        
        dif_value = sum(abs(W - W_new))/sum(abs(W))
        W = W_new
        t.0 = t_new
        count = count + 1
        res = tr( t(W) %*% Sigma %*% W) * 0.5 - tr(t(W) %*%U) + lambda *  sum (abs(W))
        if(sum(abs(W)) == 0){
            break
        }
        cat(sep = "", "[",count,"]","dif=", round(dif_value, 3), "; Loss=", res, "\n")
        
    }
    return(W)
}

CCA_GEV =  function(X, Y, lambdax = 0.1, lambday = 0.01, diff_thre = 2e-6, max_iter = 500,
                    k = 2, standardize = FALSE){
    #####################################################
    ##### Sparse CCA using GEV framework                                        #####
    ##### X is the n x p matrix                                                              #####
    ##### Y is a n x m matrix                                                                    #####
    ##### lambdax and lambday are the tuning parameters            #######
    ##### standardize indictas whether the data need to be           #######
     ##### standardized to mean 0 and sd of 1                                #######
     ##### k is the intrinsic dim
    ######################################################    
    if(standardize){
        X = Col_standardize(center_rowmeans(X))
        Y = Col_standardize(center_rowmeans(Y))
    }
    
    p = ncol(X); m = ncol(Y); n = nrow(X)
    lambda = c(rep(lambdax, p), rep(lambday, m))
    Sig = Ome = matrix(0, p + m, p + m)
    Ome[1:p, (p + 1): (p + m)] = t(X) %*% Y / n
    Ome[ (p + 1): (p + m), 1:p] = t(Ome[1:p, (p + 1): (p + m)])
    
    Sig[1:p, 1:p] = t(X) %*% X /n
    Sig[(p + 1): (m + p), (p + 1): (m + p)] = t(Y) %*% Y /n
    U_est = Get_U(Ome, k)
    W2 = GEV_FISTA(Sig, U_est, lambda = lambda, diff_thre = diff_thre, max_iter = max_iter)
    return(W2)    
}

FDA_GEV =  function(X, Y, lambda = 0.1, diff_thre = 2e-6, max_iter = 500,
        k = 1, standardize = TRUE){
    #####################################################
    ##### Sparse FDA using GEV framework                                        #####
    ##### X is the n x p matrix                                                                 ####
    ##### Y is a n x 1 label vector                                                           ####
    ##### lambda and lambday are the tuning parameters                   ####
    ##### standardize indictas whether the data need to be                 ####
    ##### standardized to mean 0 and sd of 1                                       ####
    ##### k is the intrinsic dim
    ######################################################    
    if(standardize){
        X = center_rowmeans(X)
    }    
    p = ncol(X);  n = nrow(X); K = length(unique(Y))
    lambda = rep(lambda, p)
    Sig = matrix(0, p, p)
   #Ome = matrix(0, p, p)
    for(i in 1:K){
        id = Y == (i - 1)
        sample_size = sum(id)
        mu_temp = apply(X[id, ], 2, mean)
        Xtemp = X[id, ] - rep(1, sample_size) %*% t(mu_temp)
        Sig = Sig +  t( Xtemp ) %*% Xtemp  
        #Ome = Ome + (sample_size/n) * mu_temp %*% t(mu_temp)
    }
            
    U_est = NULL
    Sig = Sig / (n - K)
    if(K == 2){
        id = Y == 1
        U_est = apply(X[id, ], 2, mean)
        a = t(U_est) %*% U_est
        U_est =  U_est / as.numeric(sqrt(a))
    }else{
        id = Y == 0
        m1 = apply(X[id, ], 2, mean)
        for(j in 1:(K-1)){
            id = Y == j
            mu_temp = apply(X[id, ], 2, mean)            
            U_est = cbind(U_est, mu_temp - m1)
        }
        U_est = U_est %*% diag(1/sqrt(diag(t(U_est) %*% U_est))   )
    }
    
    #U_est = Get_U(Ome, k)
    W = GEV_FISTA(Sig, as.matrix(U_est), lambda = lambda, diff_thre = diff_thre, max_iter = max_iter)
    return(W)    
}


FDA_pred = function(X, Y, X_test, Y_test, lambda = 0.1, diff_thre = 2e-6, max_iter = 500,
        k = 1, standardize = TRUE){
    #####################################################
    ##### Prediction Sparse FDA using GEV framework                         ####
    ##### X is the n x p matrix                                                                 ####
    ##### Y is a n x 1 label vector                                                           ####
    ##### X_test is the m x p matrix                                                          ####
    ##### Y is a m x 1 label vector                                                           ####
    ##### lambda and lambday are the tuning parameters                   ####
    ##### standardize indictas whether the data need to be                 ####
    ##### standardized to mean 0 and sd of 1                                       ####
    ##### k is the intrinsic dim                                                                  ####
    ######################################################    
    W  = FDA_GEV(X, Y, lambda = lambda, diff_thre = diff_thre,  max_iter = max_iter,
            k = k, standardize = standardize)
    Wdim = apply(abs(W), 2, sum) != 0
    if(sum(Wdim) == 0){
        result = list(error = 0.5 * length(Y_test), Ypred = 0, W = W)  
    }else{
        W = as.matrix(W[, Wdim]) # prevent column of all 0
        mtotal = apply(X, 2, mean)
        p = ncol(X);  n = nrow(X); K = length(unique(Y))
        mu = NULL
        for(i in 1:K){
            id = Y == (i - 1)
            sample_size = sum(id)
            mu = cbind(mu, apply(X[id, ], 2, mean))
        }
        Z = ( (X - rep(1, length(Y)) %*% t(mtotal)) %*% W)
        lda_res = lda(x = Z, grouping = Y)
        Z_test = ( (X_test - rep(1, length(Y_test)) %*% t(mtotal)) %*% W)     
        Ypred = predict(lda_res, Z_test)$class
        
        error = sum(Ypred != Y_test)
        result = list(error = error, Ypred = Ypred, W = W)        
    }
    return(result)
}


Orc_pred = function(mu, Sigma, X_test, Y_test){
    #####################################################
    ##### Oracle pediction Sparse FDA using GEV framework                 ####
    ##### mu is the p x K matrix containing centroid                                ####
    ##### Sigma is a p x p Covariance  matrix                                          ####
    ##### X_test is the m x p matrix                                                          ####
    ##### Y is a m x 1 label vector                                                           ####
    ######################################################        
    K = dim(mu)[2]; p = dim(mu)[1]
    U = matrix(0, p, K - 1); n = dim(X_test)[1]
    mtotal = apply(mu, 1, mean)
    Post = matrix(0, n, K)
    W = matrix(0, p, K)
   
    for(i in 1:K){
        W[, i] =  solve(Sigma) %*% (mu[, i])
        Post[, i] = (X_test - 0.5 * rep(1, n) %*% t(mu[, i])) %*% W[, i]        
        #U[, i] = solve(Sigma) %*% (mu[, i + 1] - mu[, 1])
    }
    
    max_post = apply(Post, 1, max)
   Ypred = Y_test
   for(i in 1:n){
       id = Post[i, ] == max_post[i]
       Ypred[i] = (1:K)[id]
   }
   
   result = list(Ypred = Ypred, error = sum(Ypred != (Y_test + 1) ), W = W)
   return(result)
}




Shrink 	= function(x, lambda){
    #################################
    # This is a function for soft-thresholding ###
    # x is a scalar                                            ###
    # lambda is the tuning parameter            ###
    ################################
    if(abs(x) >= lambda){
        result = sign(x) *(abs(x) - lambda)
    } else{
        result = 0
    }
    return(result)
}

center_rowmeans <- function(X) {
    # X is a n x p matrix
    # Center the variables
    xcenter = apply(X, 2, mean)
    X = X - rep(xcenter, rep.int(nrow(X), ncol(X)))
    return(X)
}


Get.projection = function(U){
    ###########################
    ### U is a p x K matrix ###
    ### Get the proj        ###
    ###########################
    Proj =  U %*% solve(t(U) %*% U) %*% t(U)
    return(Proj)
}

Get.error = function(b, b.h, X){
    #####################################
    ##### b is the true basis matrix ####
    ##### b.h is the estimate        ####
    #####################################
    P = Get.projection(b) 
    d = dim(b)[2]
    d.hat = dim(as.matrix(b.h))[2]
    B =	t(b.h) %*% b.h
    B.eigen = eigen(B)
    if(B.eigen$values[d.hat] < 0.0001)
    {
        id = B.eigen$values == 0
        B.eigen$values[id] = 0.0001
#        B.eigen$values[d.hat] = 0.0001
        B = B.eigen$vectors %*% diag(B.eigen$values) %*% t(B.eigen$vectors)
    }
    P.h = b.h%*% solve(B) %*% t(b.h)
    temp = P - P.h
    Forb = sum(diag(temp %*% t(temp)))
    Trace = tr(P %*% P.h)/ d
    Cancor = cancor(X %*% b, X %*% b.h)$cor[1]
    result = list(Forb = Forb, Trace = Trace, Cancor = Cancor)
    return(result)
}	


Get_errorF = function(b, b.h){
    #####################################
    ##### b is the true basis matrix ####
    ##### b.h is the estimate        ####
    #####################################
    P = Get.projection(b) 
    d = dim(b)[2]
    d.hat = dim(as.matrix(b.h))[2]
    B =	t(b.h) %*% b.h
    B.eigen = eigen(B)
    if(B.eigen$values[d.hat] < 0.0001)
    {
        id = B.eigen$values == 0
        B.eigen$values[id] = 0.0001
        B = B.eigen$vectors %*% diag(B.eigen$values) %*% t(B.eigen$vectors)
    }
    P.h = b.h%*% solve(B) %*% t(b.h)
    temp = P - P.h
    Forb = sum(diag(temp %*% t(temp)))^0.5
    return(Forb)
}	

cv.part <- function(n, k) {
    ntest <- floor(n/k)
    ntrain <- n-ntest
    
    ind <- sample(n)
    
    trainMat <- matrix(NA, nrow=ntrain, ncol=k)
    testMat <- matrix(NA, nrow=ntest, ncol=k)
    
    nn <- 1:n
    
    for (j in 1:k) {
        sel <- ((j-1)*ntest+1):(j*ntest)
        testMat[,j] <- ind[sel ]
        sel2 <-nn[ !(nn %in% sel) ]
        trainMat[,j] <- ind[sel2]
    }
    
    return(list(trainMat=trainMat, testMat=testMat))
}

Col_standardize = function(X){
    # standardize a matrix to have col standard dev of 1
    # X: n x p matrix
    X_sd = apply(X, 2, sd)
    p = dim(X)[2]
    for(i in 1:p){
        X[, i] = X[, i]/X_sd[i]
    }
    return(X)
}

CCACV = function (X, Y, fold = 5, lambdax= seq(0.001, 0.01, length =10), k = 2,
       lambday =  seq(0.001, 0.01, length =10), standardize = FALSE, seed = 1, diff_thre = 1e-5,
       max_iter = 1000) 
{
    ####################################
    ## CV for selecting tuning parameter lambda for CCA
    ## X: n x p matrix
    ## Y: n x m matrix 
    ## fold: fold of CV
    ## lambdax: the seq of lambda for X
    ## lambday: the seq of lambda for Y
    if(standardize){
        X = Col_standardize(center_rowmeans(X))
        Y = Col_standardize(center_rowmeans(Y))
    }
    set.seed(seed)
    
    p <- ncol(X)
    n <- nrow(X)
    m = ncol(Y)
    part.list <- cv.part(n, fold)
    lambdax_n = length(lambdax)
    lambday_n = length(lambday)
    cca.re <- array(0, c(fold,lambdax_n, lambday_n))
    for (f in 1:fold) {                    
        train.id = part.list$trainMat[, f]
        train.n  = length(train.id)
        test.id  = part.list$testMat[, f]   
        test.n   = length(test.id)
        trainX  = X[train.id, ]
        testX   = X[test.id, ]
        trainY  = Y[train.id, ]
        testY   = Y[test.id, ]  
        
        for( i in 1:lambdax_n){
            for ( j in 1:lambday_n){
                 W_train = CCA_GEV(trainX, trainY, lambdax = lambdax[i], lambday = lambday[j], k = k, 
                         diff_thre = diff_thre, max_iter = max_iter)
                 U_train = W_train[1:p, ]
                 V_train = W_train[(p+1):(p + m), ]
                 ### Prevent some direction has all zeros
                 if(min(apply( abs( testX %*% U_train), 2, sum), apply(abs(testY %*% V_train), 2, sum)) ==0){
                     cca.re[f, i, j] = 0
                 }else{
                     if(k == 1){
                         cca.re[f, i, j] =  cor(testX %*% U_train,  testY %*% V_train)
                     }else{
                     cca.re[f, i, j] =  CCA(x = testX %*% U_train, z = testY %*% V_train,  penaltyx = 0,
                             standardize = FALSE,  penaltyz = 0, K = 1)$cors
                 }
                 }
                cat(sep = "", "[fold=",f,"]","i=",i, "j=", j,"\n")
            }
        }
    }
    
    cca.mean <- apply(cca.re, c(2, 3), mean)
    cca.sd <- apply(cca.re, c(2, 3), sd)
    
    id = which(cca.mean == max(cca.mean))  # finding the largest cca
    i = (id - 1) %%lambdax_n + 1
    j = floor(( id - 1) / lambdax_n) + 1
    
    lambdaopt = c(lambdax[i] , lambday[j])
    
    outlist <- list(cca.mean = cca.mean, cca.sd = cca.sd, lambdaopt = lambdaopt)
    return(outlist)
}                


FDACV = function (X, Y, fold = 5, lambda= seq(0.001, 0.01, length =10), k = 2,
        standardize = FALSE, seed = 1, diff_thre = 1e-5, max_iter = 1000) 
{
    ####################################
    ## CV for selecting tuning parameter lambda for CCA
    ## X: n x p matrix
    ## Y: n x 1 vector  taking values of "0, 1, ..., K"
    ## fold: fold of CV
    ## lambda: the seq of lambda
    ####################################
    if(standardize){
        X = center_rowmeans(X)
    }
    set.seed(seed)
    
    p = ncol(X);  n = nrow(X); K = length(unique(Y)); 
    sample.size = rep(0, K)
    trainID = NULL
    testID = NULL
    temp_n = 0
    for(i in 1:K){
        sample.size[i] = sum(Y == (i -1))
        temp = cv.part(sample.size[i], fold)
        trainID = rbind(trainID, temp$trainMat + temp_n)
        testID = rbind(testID, temp$testMat + temp_n)
        temp_n = temp_n + sample.size[i]
    }
    
    lambda_n = length(lambda)
    err.re <- matrix(0, fold, lambda_n)
    for (f in 1:fold) {                    
        trainX  = X[trainID[, f], ]; testX   = X[testID[, f], ]
        trainY  = Y[trainID[, f] ];  testY   = Y[testID[, f] ]          
        for( i in 1:lambda_n){
                #W_train = FDA_GEV(trainX, trainY, lambda = lambda[i], k = k, diff_thre = diff_thre,
                #              max_iter = max_iter)
                #Xmean = apply(trainX, 2, mean)
                #Ypredict =( (testX - rep(1, dim(testX)[1]) %*% t(Xmean)) %*% W_train) > 0
               # err.re[f, i] =  min( sum(Ypredict == testY), sum(!Ypredict == testY)   )
                err.re[f, i] = FDA_pred(trainX, trainY, testX, testY, lambda = lambda[i],
                                    diff_thre = diff_thre, max_iter = max_iter, k = k )$error  
                cat(sep = "", "[fold=",f,"]","i=",i,  "\n")
            }
        }
        
    err.mean <- apply(err.re, 2, mean)
    err.sd <- apply(err.re, 2, sd)
    
    id = which(err.mean == min(err.mean))  # finding the largest cca
    lambdaopt = lambda[id]
    
    outlist <- list(err.mean = err.mean, err.sd = err.sd, lambdaopt = lambdaopt, fold = fold,
            lambda_seq = lambda)
    return(outlist)
}                


Get_error1 = function(a, esta){
    # estimate the sin of the two vectors a and esta
    scale_a =  as.numeric(t(a) %*% a)
    A = a %*% t(a)/ scale_a
    scale_aest =  as.numeric( t(esta) %*% esta)
    Aest = esta %*% t(esta) / scale_aest
     res = tr( (A - Aest) %*% t(A - Aest))
    return(res)
}

Get_errorP = function(A, Aest, Sig){
    # get the prediction loss
    r = ncol(A)
    for(i in 1:r){
        A[, i] = A[, i] / as.numeric(t(A[, i]) %*% Sig %*% A[, i]  )
        Aest[, i] =  sign(cor(A[, i], Aest[, i])) * Aest[, i] / as.numeric(sqrt( t(Aest[, i]) %*% Sig %*% Aest[, i]  ))
    }
    
    res = tr( t(Aest - A) %*% Sig %*%  (Aest - A))
    return(res)
}


NearestNeighborSigma1 = function(p, m,sd){
#Simulate nearest neighbor network sollow Li&Gui 2006
#p is the dimension
#m control the sparsity
#sd is seed number
    
    set.seed(sd)
    x = runif(p, 0, 1)
    y = runif(p, 0, 1)
    
    d = matrix(NA, p, p)
    Omega = matrix(0, p, p)
    
    for(i in 1:p){
        for(j in 1:p){
            d[i,j] = sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2)
        }
    }
    
    m.set = matrix(NA,p,m)
    for(i in 1:p){
        m.set[i,] = order(d[i,])[2:(m+1)]
    }
    
    for(i in 1:p){
        id = m.set[i,]
        for(j in 1:m){
            while({sum(m.set[id[j],] == i) ==1} & {Omega[i,id[j]] == 0 }){
                Omega[i,id[j]] = runif(1, 1, 1)*(-1)^rbinom(1,1,0.5)
                Omega[id[j], i] = Omega[i, id[j]]
            }
        }
        Omega[i, i] = sum(abs(Omega[i,]))
        if(Omega[i,i] == 0){Omega[i,i] = 1}
        
    }
    
    Omega.diag = diag(Omega)
    for(i in 1:p){
        Omega[i,] = Omega[i,]/sqrt(Omega.diag[i])
        Omega[,i] = Omega[,i]/sqrt(Omega.diag[i])
        
    }    
    
    result = list (Omega = Omega)
    result
}


NearestNeighborSigma2 = function(p, m, sd, min_val = 0.2, max_val = 0.8 ){
# Simulate nearest neighbor network sollow Li&Gui 2006
# p is the dimension
# m control the sparsity
# sd is seed number    
    set.seed(sd)
    x = runif(p, 0, 1)
    y = runif(p, 0, 1)
    
    d = matrix(NA, p, p)
    Omega = matrix(0, p, p)
    
    for(i in 1:p){
        for(j in 1:p){
            d[i, j] = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
        }
    }
    
    m.set = matrix(NA, p, m)
    for(i in 1:p){
        m.set[i,] = order(d[i, ])[2:(m + 1)]
    }
    
    for(i in 1:p){
        id = m.set[i,]
        for(j in 1:m){
            while({sum(m.set[id[j],] == i) ==1} & {Omega[i,id[j]] == 0 }){
                Omega[i,id[j]] = runif(1, min_val, max_val)*(-1)^rbinom(1, 1, 0.5)
                Omega[id[j], i] = Omega[i, id[j]]
            }
        }
        Omega[i, i] = sum(abs(Omega[i,]))
        if(Omega[i,i] == 0){Omega[i,i] = 1}
    }
    Omega.diag = diag(Omega)
    for(i in 1:p){
        Omega[i,] = Omega[i,]/sqrt(Omega.diag[i])
        Omega[,i] = Omega[,i]/sqrt(Omega.diag[i])        
    }    
    Sig = solve(Omega)
    Sigma = diag(1/ sqrt(diag(Sig))) %*% Sig %*% diag(1/ sqrt(diag(Sig))) 
    result = list (Omega = Omega, Sigma = Sigma)
    return(result)
}


