dsda.path<-function(x.matrix,y,standardize=FALSE,lambda=lambda,alpha=1,thresh=1e-7){
     n<-length(y)
     n1<-n-sum(y)
     n2<-sum(y)
 
     obj<-glmnet(x.matrix,y,standardize=standardize,family="gaussian",alpha=alpha,thresh=thresh)
     if(missing(lambda))lambda<-obj$lambda
     beta<-coef(obj,s=lambda,standardize=standardize)
      beta1<-as.matrix(beta[-1,,drop=FALSE])
      sel<-apply(as.matrix(abs(beta1)),1,sum)
      sel<-which(sel!=0)
      mu1<-as.matrix(apply(x.matrix[y==0,sel,drop=F],2,mean))
      mu2<-as.matrix(apply(x.matrix[y==1,sel,drop=F],2,mean))
      sigma.hat<-(cov(x.matrix[y==0,sel,drop=F])*(n1-1)+cov(x.matrix[y==1,sel,drop=F])*(n2-1))/(n-2)
     beta1[sel,]<-sweep(as.matrix(beta1[sel,,drop=F]),2,sign(t(mu2-mu1)%*%beta1[sel,,drop=F]),"*")
     beta0<--t((mu1+mu2))%*%beta1[sel,,drop=F]/2-log(n1/n2)*diag(t(beta1[sel,,drop=F])%*%sigma.hat%*%beta1[sel,,drop=F])/(t(mu2-mu1)%*%beta1[sel,,drop=F])
     beta0[is.na(beta0)]<-n2-n1
     beta[1,]<-beta0

     list(beta=beta,lambda=lambda)
}
