\name{dsda.path}
\alias{dsda.path}
\title{Solution path for direct sparse discriminant analysis}
\description{
Compute the solution path for direct sparse discriminant analysis.}
\usage{
dsda.path(x.matrix,y,standardize=FALSE,lambda=lambda,alpha=1,thresh=1e-7)
}
\arguments{
\item{x.matrix}{an n by p matrix containing the predictors.}
\item{y}{an n-dimensional vector containing the class labels. The classes have to be labeled as 0 and 1.}
\item{standardize}{a logic object indicating whether x.matrix should be standardized before performing DSDA. Default is FALSE.}
\item{lambda}{a sequence of lambda's.}
\item{alpha}{the elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{thresh}{convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}
\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}
\value{
\item{beta}{the solution path.}
\item{lambda}{the sequence of lambda's used in computing the solution path.}
}
\examples{
    data("prostate")
    x<-prostate.x
    y<-prostate.y
    obj.path <- dsda.path(x, y)

}
