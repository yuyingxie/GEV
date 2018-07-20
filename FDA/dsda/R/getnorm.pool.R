getnorm.pool<-function(x,y){
  f0<-ecdf(x[y==0])
  f1<-ecdf(x[y==1])
  n0<-sum(y==0)
  n1<-sum(y==1)
  n<-n0+n1
  delta.n0<-1/n0^2
  delta.n1<-1/n1^2
  
  v0<-rep(0,n)
  v0[y==0]<-qnorm(truncate(x[y==0],f0,delta.n0))
  v0[y==1]<-qnorm(truncate(x[y==1],f0,delta.n0))
  mu0.hat<-mean(v0[y==1])
  
  v1<-rep(0,n)
  v1[y==0]<-qnorm(truncate(x[y==0],f1,delta.n1))
  v1[y==1]<-qnorm(truncate(x[y==1],f1,delta.n1))
  mu1.hat<-mean(v1[y==0])
  mu.hat<-n0/n*mu0.hat-n1/n*mu1.hat  
  
  transform<-function(t){
    n0*qnorm(truncate(t,f0,delta.n0))/n+n1*(qnorm(truncate(t,f1,delta.n1))+mu.hat)/n
  }
  x.norm<-transform(x)
  
  list(x.norm=x.norm,f0=f0,f1=f1,mu.hat=mu.hat,transform=transform)
}