truncate<-function(x,f,delta.n){
  y<-f(x)
  y[y<delta.n]<-delta.n
  y[y>1-delta.n]<-1-delta.n
  y
}
