getnorm<-function(x,y,type="pooled"){
 type<-match.arg(type,c("pooled","naive"))
 if(type=="pooled"){getnorm.pool(x,y)}else{
      getnorm.naive(x,y)
  }
}

