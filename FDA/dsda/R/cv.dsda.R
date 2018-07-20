 cv.dsda<-function (x, y, K = 10, standardize=FALSE,alpha=1,lambda=lambda,thresh=1e-7,lambda.opt="min")
 {
     n<-length(y)
     n1<-n-sum(y)
     n2<-sum(y)

     if(K>n)stop("The number of folds should be smaller than the sample size.")

     all.folds <- cv.folds(length(y), K)
     
     if(missing(lambda)){
        fit <- glmnet(x, y,  family="gaussian",alpha=alpha,standardize=standardize,thresh=thresh)
        lambda<-fit$lambda}

     nlambda<-length(lambda)
     residmat <- matrix(0, nlambda, K)
     for (i in seq(K)) {
         omit <- all.folds[[i]]
         fit <- dsda.path(x[-omit, , drop = FALSE], y[-omit], lambda=lambda,standardize=standardize,thresh=thresh)
         fit<-predict.dsda(fit,x[omit,,drop=FALSE])
         residmat[, i] <- apply(abs(sweep(fit,1,y[omit])),2,mean)}
       
     residmat[is.na(residmat)]<-min(n1/n,n2/n)
     residmat<-matrix(residmat,nrow=nlambda)
     cv <- apply(residmat, 1, mean)
     cv.error <- sqrt(apply(residmat, 1, var)/K)
     if(lambda.opt=="min"){
        bestlambda<-max(lambda[which(cv==min(cv))])}
     else{
        bestlambda<-max(lambda[which(cv==min(cv))])}
     object <- list(lambda = lambda, cv = cv, cv.error = cv.error,bestlambda=bestlambda)
     invisible(object)
 }
 
