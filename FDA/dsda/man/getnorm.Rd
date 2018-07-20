\name{getnorm}
\alias{getnorm}
\title{Direct sparse discriminant analysis}
\description{
Transform the predictors to achieve normality.
}
\usage{
getnorm(x,y)
}
\arguments{
\item{x}{an n dimensional vector containing n observations for one predictor.}
\item{y}{an n-dimensional vector containing the class labels.}
\item{type}{The type of estimator. Two estimators were proposed in Mai & Zou (2015), the naive estimator and the pooled estimator. The function getnorm() uses the naive estimator if type="naive", and it uses the pooled estimator if type="pool". The default is "pool". When the naive estimator is used, it is recommended to label the class with more samples as Class 0.}
}
\value{
\item{x.norm}{Transformed x.}
\item{f0}{The transformation computed based on observations from Class 0. Not applicable if type="naive".}
\item{f1}{The transformation computed based on observations from Class 1. Not applicable if type="naive".}
\item{mu.hat}{The sample mean for transformed x from Class 1.}
\item{transform}{The transformation that was actually used to transform x.}
}
\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.

Mai, Q. and Zou, H. (2015). Sparse semiparametric discriminant analysis. Journal of Multivariate Analysis, 135, 175-188. 
}

\examples{
data("malaria")   ##load the malaria data
x<-malaria.x
y<-malaria.y


n<-length(y)
n1<-sum(y)
n2<-n-n1
n1.test<-round(n1/2)
n2.test<-round(n2/2)
n.test<-n1.test+n2.test
n.train<-n-n.test
id.test<-c(sample(which(y==1),n1.test),sample(which(y==0),n2.test))

p<-ncol(x)
x.train<-x[-id.test,]
y.train<-y[-id.test]
x.test<-x[id.test,]
y.test<-y[id.test]

 x.norm<-matrix(0,n.train,p)
  x.test.norm<-matrix(0,n.test,p)
  for(i in 1:p){
    obj.norm<-getnorm(x.train[,i],y.train)
    x.norm[,i]<-obj.norm$x.norm
    x.test.norm[,i]<-obj.norm$transform(x.test[,i])
  }
    
  obj<-dsda(x.norm,y.train,x.test.norm,y.test)



}