predict.dsda<-function(obj,x.matrix){
  beta<-obj$beta
  pred<-as.matrix(cbind(1,x.matrix)%*%beta)
  pred<-ifelse(pred>0,1,0)
  invisible(pred)
}
