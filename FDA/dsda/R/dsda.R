dsda<-function(x.matrix,y,x.test.matrix=NULL,y.test=NULL,standardize=FALSE,lambda.opt="min",K=10,lambda=lambda,alpha=1,thresh=1e-7){
 n<-nrow(x.matrix)
 d<-ncol(x.matrix)
 n2<-sum(y)
 n1<-n-n2
 
 la<-cv.dsda(x.matrix,y,standardize=standardize,alpha=alpha,K=K,lambda=lambda,thresh=thresh,lambda.opt=lambda.opt) 
 s<-la$bestlambda
 
 obj.path<-dsda.path(x.matrix,y,lambda=s,standardize=standardize,thresh=thresh)
 beta<-obj.path$beta
 
 if(!missing(x.test.matrix)){
   n.test<-dim(x.test.matrix)[1]
   pred<-predict.dsda(obj.path,x.test.matrix)
   error<-mean(pred!=y.test)
   }else{error<-NA}
 result<-list(error=error,beta=beta,s=s)}
