getnorm.naive<-function(x,y){
  f0<-ecdf(x[y==0])
  n0<-sum(y==0)
  n1<-sum(y==1)
  n<-n0+n1
  delta.n0<-1/n0^2
  
  v0<-rep(0,n)
  v0[y==0]<-qnorm(truncate(x[y==0],f0,delta.n0))
  v0[y==1]<-qnorm(truncate(x[y==1],f0,delta.n0))
  mu0.hat<-mean(v0[y==1])
  
  
  transform<-function(t){
    qnorm(truncate(t,f0,delta.n0))  }
  x.norm<-transform(x)
  
  list(x.norm=x.norm,mu.hat=mu0.hat,transform=transform)
}