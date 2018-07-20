\name{cv.dsda}
\alias{cv.dsda}
\title{Cross validation for direct sparse discriminant analysis}
\description{
Choose the optimal lambda for direct sparse discriminant analysis by cross validation.
}
\usage{
cv.dsda(x, y, K = 10, standardize=FALSE,alpha=1,lambda=lambda,thresh=1e-7,lambda.opt="min")
}
\arguments{
\item{x}{an n by p matrix containing the predictors.}
\item{y}{an n-dimensional vector containing the class labels.}
\item{standardize}{a logic object indicating whether x.matrix should be standardized before performing DSDA. Default is FALSE.}
\item{lambda.opt}{should be either "min" or "max", specifying whether the smallest or the largest lambda with the smallest cross validation error should be used for the final classification rule.}
\item{K}{the number of folds to be used in cross validation. Default is 10.}
\item{lambda}{a sequence of lambda's.}
\item{alpha}{the elasticnet mixing parameter, the same as in glmnet. Default is alpha=1 so that the lasso penalty is used.}
\item{thresh}{convergence threshold for coordinate descent, the same as in glmnet. Default is 1e-7.}
}

\references{
Mai, Q., Zou, H. and Yuan, M. (2013). A direct approach to sparse discriminant analysis in ultra-high dimensions. Biometrika, 99, 29-42.
}

\seealso{
\code{cv.dsda}
\code{predict.dsda}
\code{dsda.path}
}
\value{
\item{lambda}{the sequence of lambda's used in cross validation.}
\item{cv}{cross validation errors.}
\item{cv.error}{the standard error of the cross validation errors.}
\item{bestlambda}{the optimal lambda chosen by cross validation.}
}
\examples{
    data("prostate")
    x<-prostate.x
    y<-prostate.y
    obj.cv <- cv.dsda(x, y)
    s <- obj.cv$bestlambda
    obj.path <- dsda.path(x, y, lambda = s)

}