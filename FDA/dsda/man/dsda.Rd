\name{dsda}
\alias{dsda}
\title{Direct sparse discriminant analysis}
\description{
Perform direct sparse discriminant analysis, with the optimal lambda chosen by cross validation.
}
\usage{
dsda(x.matrix,y,x.test.matrix=NULL,y.test=NULL,standardize=FALSE,lambda.opt="min",K=10,lambda=lambda,alpha=1,thresh=1e-7)
}
\arguments{
\item{x.matrix}{an n by p matrix containing the predictors.}
\item{y}{an n-dimensional vector containing the class labels.}
\item{x.test.matrix}{(optional.) the predictors of a testing set.}
\item{y.test}{(required if x.test.matrix is supplied, but otherwise optional.) the class labels of the testing set.}
\item{standardize}{a logic object indicating whether x.matrix should be standardized before performing DSDA. Default is FALSE.}
\item{lambda.opt}{should be either "min" or "max", specifying whether the smallest or the largest lambda with the smallest cross validation error should be used for the final classification rule.}
\item{K}{the number of folds to be used in cross validation. Default is 10.}
\item{lambda}{a sequence of lambda's.}
\item{alpha}{the elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{thresh}{convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}

\references{
A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}
\seealso{
\code{cv.dsda}
\code{predict.dsda}
\code{dsda.path}
}
\value{
\item{error}{testing error if x.test.matrix is supplied.}
\item{beta}{the coefficients of the classification rule corresponding to the optimal lambda chosen by cross validation.}
\item{s}{the optimal lambda chosen by cross validation.}
}
\examples{
data("prostate")   ##load the prostate data
x<-prostate.x
y<-prostate.y

n<-length(y)            ##split the original dataset to a training set and a testing set
n.test<-round(n/3)
set.seed(20120822)
id<-sample(n,n.test,replace=FALSE)
x.train<-x[-id,]
x.test<-x[id,]
y.train<-y[-id]
y.test<-y[id]

set.seed(123)
obj<-dsda(x.train,y.train,x.test,y.test)  ##perform direct sparse discriminant analysis
obj$error 
}